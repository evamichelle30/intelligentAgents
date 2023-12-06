package template;

import logist.LogistPlatform;
import logist.LogistSettings;
import logist.agent.Agent;
import logist.behavior.AuctionBehavior;
import logist.config.Parsers;
import logist.plan.Plan;
import logist.simulation.Vehicle;
import logist.task.DefaultTaskDistribution;
import logist.task.Task;
import logist.task.TaskDistribution;
import logist.task.TaskSet;
import logist.topology.Topology;
import logist.topology.Topology.City;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


@SuppressWarnings("unused")
public class AuctionTemplate implements AuctionBehavior {

	private Topology topology;
	private DefaultTaskDistribution distribution;
	private Agent agent;
	private double p; // probability of returning old solution for SLS algorithm
	private long timeout_setup;
	private long timeout_bid;
	private long timeout_plan;
	private Random random;

	private Candidate mySolution;
	private Candidate opponentSolution;

	private double opponentBidEstimate;
	private double maxMarginalCost;
	private double maxCapacity;
	private double opponentCostMultiplier = 1;

	@Override
	public void setup(Topology topology, TaskDistribution distribution,
			Agent agent) {

		this.topology = topology;
		this.distribution = (DefaultTaskDistribution) distribution;
		this.agent = agent;
		this.maxMarginalCost = 0;

		this.maxCapacity = 0;
		for (Vehicle vehicle : agent.vehicles()) {
			this.maxCapacity = Math.max(vehicle.capacity(), maxCapacity);
		}

		long seed = -9019554669489983951L * agent.id();
		this.random = new Random(seed);

		this.p = 0.2; // set p

		this.mySolution = new Candidate(agent.vehicles());
        this.opponentSolution = new Candidate(agent.vehicles());

		// this code is used to get the timeouts
        LogistSettings ls = null;
        try {
            ls = Parsers.parseSettings("config" + File.separator + "settings_auction.xml");
        } catch (Exception exc) {
            System.out.println("There was a problem loading the configuration file.");
        }

		long timeout_margin = 200;

		// the setup method cannot last more than timeout_setup milliseconds
		timeout_setup = LogistPlatform.getSettings().get(LogistSettings.TimeoutKey.SETUP) - timeout_margin;
		// the bid method cannot last more than timeout_bid milliseconds
		timeout_bid = LogistPlatform.getSettings().get(LogistSettings.TimeoutKey.BID) - timeout_margin;
		// the plan method cannot last more than timeout_plan milliseconds
		timeout_plan = LogistPlatform.getSettings().get(LogistSettings.TimeoutKey.PLAN) - timeout_margin;
		System.out.println("Agent " + agent.id() + ": timeout_setup: " + timeout_setup + ", timeout_bid: " + timeout_bid + ", timeout_plan: " + timeout_plan);
	}

	@Override
	public void auctionResult(Task previous, int winner, Long[] bids) {
		if (winner == agent.id()) {
            System.out.println("You won task: " + previous.id);
            mySolution.addTask(previous);
        } else {
        	System.out.println("Opponent won task: " + previous.id);
            opponentSolution.addTask(previous);
        }

        // Calibrate my estimate of the opponent's bid (assuming we have only 1 opponent)
		long opponentBid = 0;
        for (int i = 0; i < bids.length; i++) {
            if (i != agent.id()) {
                opponentBid = bids[i];
                updateOpponentCostMultiplier(opponentBid);
            }
        }
	}


	public void updateOpponentCostMultiplier(long opponentBid) {
		double discount = 0.5;
		if (opponentBidEstimate == 0)
			return;
		double ratio = opponentBid / opponentBidEstimate;
		opponentCostMultiplier = discount * ratio + (1 - discount) * opponentCostMultiplier;
		return;
	}

	
	/**
	 * Naive bidding strategy that takes into account (i) the marginal cost, (ii) a crude modeling of the opponent, and (iii) potential future tasks.
	 * 
	 * To compute the marginal cost, we use the Stochastic Local Search (SLS) algorithm implemented in the previous exercise (see marginalCostEstimation).
	 * Of course you may substitute this with a more sophisticated optimization technique (e.g., the simplex algorithm).
	 * 
	 * To model the opponent we assume the same vehicles and starting cities, and then try to "learn" a multiplier (opponentCostMultiplier) to fit the estimated opponent's bid to the actual bid (see updateOpponentCostMultiplier).
	 * More sophisticated approaches may be used (e.g., using the history of bids to update your belief on the number of vehicles and their initial location).
	 * Alternatively you may opt to not model the opponent, and use the additional computation time to optimize your own plan and compute a more accurate marginal cost.
	 * 
	 * Finally, we use the task probability distribution to sample tasks and estimate potential cost savings due to spatial overlap (see estimatePotentialForFutureCostReduction).
	 * This is a really simple / naive approach and requires to run SLS multiple times. As an alternative you could compute the probability of a task with spatial overlap (source/destination path overlapping with your current plan) being auctioned in the future.
	 * 
	 */
	@Override
	public Long askPrice(Task task) {
		// return null if we don't have a vehicle big enough to carry the task
		if (task.weight > this.maxCapacity) {
			return null;
		}

		double upper_bidding_limit = maxMarginalCost * 10;

		double margin = 1.5;
		double discount = 0.1;
		System.out.println("Auctioning task: " + task);

		List<Task> additionalTasks = new ArrayList<Task>();
		additionalTasks.add(task);


		// TODO: time out riskier?
		double myMarginalCost = marginalCostEstimation(mySolution, additionalTasks, timeout_bid * 0.25);
		double opponentMarginalCost = marginalCostEstimation(opponentSolution, additionalTasks, timeout_bid * 0.25);

		this.maxMarginalCost = Math.max(myMarginalCost, maxMarginalCost);

		// makes us loose money...
		double potentialFutureSavings = Math.min(estimatePotentialForFutureCostReduction(mySolution, task, 1, timeout_bid * 0.25), 0);
		double myCost = Math.max(0, myMarginalCost + discount * potentialFutureSavings); // potentialFutureSavings is <= 0

		potentialFutureSavings = Math.min(estimatePotentialForFutureCostReduction(opponentSolution, task, 1, timeout_bid * 0.25), 0);
		double opponentCost = Math.max(0, opponentMarginalCost + discount * potentialFutureSavings);

		opponentCost = Math.round(Math.max(opponentCostMultiplier * opponentCost, 0));
		opponentBidEstimate = opponentCost;

		System.out.println("My cost estimate: " + myCost);
		System.out.println("Opponent's cost estimate: " + opponentCost);

		double bid;
		if (opponentCost < myCost)
			bid = myCost;
        else
        	bid = myCost + discount * (opponentCost - myCost);
        // TODO: why negative reward
        if (bid < 100)
        	bid = 250 + random.nextInt(500);

        System.out.println("My bid: " + bid);

        return Math.round(Math.max(bid, 0)*margin);
	}

	public double marginalCostEstimation(Candidate solution, List<Task> additionalTasks, double timeout) {
		Candidate augmentedSolution = new Candidate(solution);
		for (Task task : additionalTasks)
			augmentedSolution.addTask(task);
		
		augmentedSolution = runSLS(augmentedSolution, timeout);
		return Math.max(augmentedSolution.cost - solution.cost, 0);
	}

	private Candidate runSLS(Candidate solution, double timeout) {
		// Begin SLS Algorithm
		long time_start = System.currentTimeMillis();

		// initial solution
		Candidate A = new Candidate(solution);

		// Optimization loop - repeat until timeout
		boolean timeout_reached = false;

		while(!timeout_reached)	{
			// record old solution
			Candidate A_old = A;	

			// generate neighbours
			List<Candidate> N = A_old.ChooseNeighbours(random);

			// Get the solution for the next iteration
			A = LocalChoice(N, A_old);

			// Check timeout condition
			if( System.currentTimeMillis() - time_start > timeout ) {
				timeout_reached = true;
			}
		}
		// End SLS Algorithm

		return A;
	}


	// Local choice to choose the next solution from the neighbours and the current solution
    public Candidate LocalChoice(List<Candidate> N, Candidate A) {
    	
    	if (random.nextFloat() < p) {	// Return A with probability p
    		
    		return A;
    		
    	}
    	else {	// Return the best neighbour with probability 1-p
    		
    		int best_cost_index = 0; // index of the neighbour with best cost until now
    		double best_cost = N.get(best_cost_index).cost; // cost of the neighbour with best cost until now

    		
    		for (int n_ind = 1; n_ind < N.size(); n_ind++ ) {
    		
    			// check if current alternative has lower cost than the current best
    			if( N.get(n_ind).cost < best_cost )	{
    				// if so, update the best solution
    				best_cost_index = n_ind;
    				best_cost = N.get(best_cost_index).cost;
    			}
    			
    		}
    			
    		// return the best solution
    		return N.get(best_cost_index);
        }
    }

    // TODO: really?
	public double estimatePotentialForFutureCostReduction(Candidate solution, Task task, int numberOfFutureTasks, double timeout) {
		int numberOfRuns = 4;
		double savings = 0;
		Candidate augmentedSolution = new Candidate(solution);
		augmentedSolution.addTask(task);
		for (int i = 0; i < numberOfRuns; i++) {
			List<Task> additionalTasks = new ArrayList<Task>();
			for (int j = 0; j < numberOfFutureTasks; j++) {
				Task potentialTask = distribution.createTask(); // Sample potential future tasks
				additionalTasks.add(potentialTask);
			}

			double marginalCostSolution = marginalCostEstimation(solution, additionalTasks, timeout / numberOfRuns / 2); // Compute the marginal cost of including the sampled tasks in your current plan.
			double marginalCostAugmentedSolution = marginalCostEstimation(augmentedSolution, additionalTasks, timeout / numberOfRuns / 2); // Compute the marginal cost of including the sampled tasks in your plan if you win the auctioned task.

			savings += (marginalCostAugmentedSolution - marginalCostSolution) / numberOfRuns;
		}

		return savings;
	}


	@Override
	public List<Plan> plan(List<Vehicle> vehicles, TaskSet tasks) {
		mySolution = runSLS(mySolution, timeout_plan);
		// Build plans for vehicles from the found solution
        //List<Plan> plans = PlanFromSolution(mySolution);
		List<Plan> plans = PlanFromSolution(mySolution);
        System.out.println("Final cost: " + mySolution.cost);
        return plans;
	}

	public List<Plan> PlanFromSolutionA(Candidate A) {
		List<Plan> planList = new ArrayList<>();

		for (int vehicleIndex = 0; vehicleIndex < A.vehicles.size(); vehicleIndex++) {
			Vehicle vehicle = A.vehicles.get(vehicleIndex);
			List<PD_Action> plan = A.plans.get(vehicleIndex);
			City currentCity = vehicle.getCurrentCity();
			Plan vehiclePlan = new Plan(currentCity);

			for (PD_Action action : plan) {
				City nextCity;
				if (action.is_pickup) {
					nextCity = action.task.pickupCity;
				} else {
					nextCity = action.task.deliveryCity;
				}

				List<City> path = AStar.findShortestPath(currentCity, nextCity);

				if (path != null) {
					for (City moveCity : path) {
						vehiclePlan.appendMove(moveCity);
					}
				}

				if (action.is_pickup) {
					vehiclePlan.appendPickup(action.task);
				} else {
					vehiclePlan.appendDelivery(action.task);
				}

				currentCity = nextCity;
			}

			planList.add(vehiclePlan);
		}

		return planList;
	}


	// Build the plan for logist platform from the candidate solution
    public List<Plan> PlanFromSolution(Candidate A) {
    	
    	// System.out.println("Constructing plan from solution...");
    	
       List<Plan> plan_list = new ArrayList<Plan>();	// create empty list of plans
      
       // Build plan for each vehicle
       for (int vehicle_ind = 0; vehicle_ind < A.vehicles.size(); vehicle_ind++) {
           
       	   Vehicle v = A.vehicles.get(vehicle_ind);
       	
       	   // get constructed plan of the vehicle
           List<PD_Action> plan = A.plans.get(vehicle_ind);
           
           // follow vehicle cities to construct plan
           City current_city = v.getCurrentCity();
           Plan v_plan = new Plan(current_city);
           
           // Append required primitive actions for each pickup/delivery action
           for (PD_Action act : plan) {
        	   
        	   City next_city;
    	   	   if(act.is_pickup) {
    	   		   next_city = act.task.pickupCity;
    	   	   }
    	   	   else {
    	   		   next_city = act.task.deliveryCity;
    	   	   }

               // Append move actions
               for(City move_city : current_city.pathTo(next_city)) {
            	   v_plan.appendMove(move_city);
               }
               // Append pickup-delivery actions
               if (act.is_pickup) {
               	v_plan.appendPickup(act.task);
               } else {
               	v_plan.appendDelivery(act.task);
               }
               current_city = next_city;
           }
           
           // add plan to the list of plans
           plan_list.add(v_plan);
       }
       return plan_list;
    }
}
