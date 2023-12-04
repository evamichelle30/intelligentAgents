package template;

//the list of imports
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import logist.LogistSettings;
import logist.Measures;
import logist.behavior.AuctionBehavior;
import logist.agent.Agent;
import logist.config.Parsers;
import logist.simulation.Vehicle;
import logist.plan.Plan;
import logist.task.Task;
import logist.task.TaskDistribution;
import logist.task.TaskSet;
import logist.topology.Topology;
import logist.topology.Topology.City;

import static algorithms.AStar.aStarPlan;

/**
 * A very simple auction agent that assigns all tasks to its first vehicle and
 * handles them sequentially.
 * 
 */
@SuppressWarnings("unused")
public class AuctionTemplate implements AuctionBehavior {

	private Topology topology;
	private TaskDistribution distribution;
	private Agent agent;
	private Random random;
	private List<Vehicle> vehicles;
	private City currentCity;
	private long currentCost;
	private List<Task> currentTasks;
	private long averageCostPerKm;
	private List<List<PD_Action>> plans;

	private long timeout_setup;
	private long timeout_plan;
	private double p;
	private double equilibriumSpread;
	private double equilibriumDiscount;

	@Override
	public void setup(Topology topology, TaskDistribution distribution,
			Agent agent) {

		this.topology = topology;
		this.distribution = distribution;
		this.agent = agent;
		this.vehicles = agent.vehicles();
		this.currentCity = vehicles.get(0).homeCity();
        this.currentTasks = (agent.getTasks() != null) ? agent.getTasks().stream().toList() : new ArrayList<>();

        long seed = -9019554669489983951L * currentCity.hashCode() * agent.id();
		this.random = new Random(seed);
		// Initialize the average cost per km.
		long totalCost = 0;
		for (Vehicle v : agent.vehicles()) {
			totalCost += v.costPerKm();
		}
		this.averageCostPerKm = totalCost / agent.vehicles().size();
		this.currentCost = 0;

		// this code is used to get the timeouts
		LogistSettings ls = null;
		try {
			ls = Parsers.parseSettings("config" + File.separator + "settings_auction.xml");
		} catch (Exception exc) {
			System.out.println("There was a problem loading the configuration file.");
		}

		// the setup method cannot last more than timeout_setup milliseconds
		timeout_setup = ls.get(LogistSettings.TimeoutKey.SETUP);

		// the plan method cannot execute more than timeout_plan milliseconds
		timeout_plan = ls.get(LogistSettings.TimeoutKey.PLAN) - 200;	// We add a little safety margin

		this.equilibriumSpread = 2.0;  // Example value, replace with your desired value
		this.equilibriumDiscount = 0.9; // Example value, replace with your desired value
	}

	@Override
	public void auctionResult(Task previous, int winner, Long[] bids) {
		if (winner == agent.id()) {
            currentTasks.add(previous);
            // update cost
			currentCost = calculatePlanCost(this.currentTasks);
		}
		// TODO: bids can contain null
		// two-parameter Weibull distribution with spread parameter m
		// TODO: update evaluation based on bids - enough?
		this.equilibriumSpread = 2.0;  // Example value, replace with your desired value
		this.equilibriumDiscount = 0.9;
	}

	@Override
	public Long askPrice(Task task) {
		// TODO: decide over all of our vehicles
		if (vehicles.get(0).capacity() < task.weight)
			return null;

		// ALWAYS BID - if not interesting bid will just be very high - else reject = null

		// Calculate equilibrium bid using the given equilibrium policy
		double equilibriumBid = calculateEquilibriumBid(task);

		return (long) Math.round(equilibriumBid);
	}

	private double calculateRationalBid(Task task) {
		long distanceTask = task.pickupCity.distanceUnitsTo(task.deliveryCity);
		long distanceSum = distanceTask
				+ currentCity.distanceUnitsTo(task.pickupCity);

		double ratio = 1.0 + (random.nextDouble() * 0.05 * task.id);
		// double bid = ratio * marginalCost;

		// Get the estimated marginal cost of the task
		double marginalCost = getMarginalCost(task);

		// Calculate the bid based on your cost
		double bid = calculateRationalBidForTask(task, (long) marginalCost);

		return (long) Math.round(bid);
	}

	private double calculateEquilibriumBid(Task task) {
		// Extract relevant parameters
		double n = agent.vehicles().size(); // TODO: update to Number of bidders
		double c = averageCostPerKm;    // Cost of doing a job

		double equilibriumBid = calculateRationalBid(task);

		// Calculate Qi and pi
		double Qi = Math.pow((n - 1) * Math.pow(equilibriumSpread, -1 / equilibriumSpread), -equilibriumSpread);
		double pi = 1 / (1 + (n - 1) * Math.pow(equilibriumBid / Qi, equilibriumSpread));

		// Calculate expected profit using the equilibrium policy formula
		double expectedProfit = c * equilibriumBid * (pi * Math.pow(equilibriumBid / Qi, -1 / equilibriumSpread) - 1);
		System.out.println("Bid: " + equilibriumBid);

		return equilibriumBid;
	}

	// TODO: difference here: we do not know tasks...
	// BUT can still be used for evaluation of a taskset!
	// Idea: before bidding - evaluate cost of new Plan (including task) - bid more than that
	// - see how much other agents got paid for actions you did not get - update evaluation
	// TODO: theoretically plan has already been computed for bidding purposes...
	@Override
	public List<Plan> plan(List<Vehicle> vehicles, TaskSet tasks) {

		System.out.println("Building plan...");

		long time_start = System.currentTimeMillis();

		// Initialize list of tasks
		List<Task> task_list = new ArrayList<>(tasks);

		// Begin SLS Algorithm - IDEA: use SLS to shift tasks and A* to do the plan per vehicles

		// create initial solution
		Candidate A = Candidate.SelectInitialSolution(random, vehicles, task_list);

		// Optimization loop - repeat until timeout TODO: can't do here, takes waaay to much time...
		boolean timeout_reached = false;

		while(!timeout_reached)	{

			// record old solution
			Candidate A_old = A;

			// generate neighbours
			List<Candidate> N = A_old.ChooseNeighbours(random);

			// Get the solution for the next iteration
			A = LocalChoice(N, A_old);

			// Check timeout condition
			if( System.currentTimeMillis() - time_start > timeout_plan ) {
				timeout_reached = true;
			}
		}

		// End SLS Algorithm


		// Build plans for vehicles from the found solution
		// TODO: use A* instead of SLS?
		List<Plan> plan = PlanFromSolution(A);

		// Informative outputs
		long time_end = System.currentTimeMillis();
		long duration = time_end - time_start;
		double cost_plan  = A.cost;

		System.out.println("The plan was generated in " + duration + " ms with a cost of " + A.cost);

		return plan;
	}

	// Local choice to choose the next solution from the neighbours and the current solution
	public Candidate LocalChoice(List<Candidate> N, Candidate A) {

		if (random.nextFloat() < p) {	// Return A with probability p
			return A;
		}
		else {	// Return the best neighbour with probability 1-p

			int best_cost_index = 0; // index of the neighbour with best cost until now
			double best_cost = N.get(best_cost_index).cost; // cost of the neighbour with best cost until now

			for (int n_ind = 1; n_ind < N.size(); n_ind++) {

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

	public List<Plan> PlanFromSolution(Candidate A) {

		System.out.println("Constructing plan from solution...");

		List<Plan> plan_list = new ArrayList<>();	// create empty list of plans

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

	private long calculateRationalBidForTask(Task task, long marginalCost) {
		// Bidding strategy based on marginal cost, speculation, and opponent modeling

		// Base bid is the marginal cost with some profit margin
		double profitMargin = 1.1; // TODO: learn margin?
		long bid = (long) (marginalCost * profitMargin);

		return bid;
	}

	private long getMarginalCost(Task newTask) {
		// Calculate cost of current plan
		long costWithoutNewTask = currentCost;

		// Create a hypothetical plan with the new task and calculate the cost
		List<Task> withNewTask =  new ArrayList<>(this.currentTasks);
        withNewTask.add(newTask);
		long costWithNewTask = calculatePlanCost(withNewTask);

		// The marginal cost of the new task
		return costWithNewTask - costWithoutNewTask;
	}

	// TODO: update - relation Candidate!
	public State generateInitialState (Vehicle vehicle, TaskSet available) {
		if (currentTasks == null || this.currentTasks.isEmpty()) {
			return new State (vehicle, available);
		} else {
			TaskSet stateCarriedTasks = (TaskSet) this.currentTasks;
			this.currentTasks.clear();
			return new State(vehicle, available, stateCarriedTasks);
		}
	}

	// TODO: maybe change to create from oldCandidate? Does saving old states bring us an advantage?
	public Candidate stateToCandidate(State state, List<Task> tasks, Vehicle vehicle) {

		int num_vehicles = vehicles.size();

		List<List<PD_Action>> plans = new ArrayList<>();
		List<List<Task>> taskLists = new ArrayList<>();
		List<Task> allTasks = new ArrayList<>(tasks);


// initialize plans and task list
		for (int i = 0; i < num_vehicles; i++) {
			plans.add(new ArrayList<>());
			taskLists.add(new ArrayList<>());
		}


// Assign all the tasks to the largest vehicle
		for (Task t : state.getCarriedTasks()) {

			List<PD_Action> plan = plans.get(vehicles.indexOf(vehicle));
			List<Task> tasks_vehicle = taskLists.get(vehicles.indexOf(vehicle));

			// Add tasks to the end of current plan
			plan.add(new PD_Action(true, t));
			plan.add(new PD_Action(false, t));

			tasks_vehicle.add(t);
		}

// calculate the cost of initial candidate solution
		double initial_cost = 0.0;
		// accumulate the cost borne by each vehicle
		for (int i = 0; i < vehicles.size(); i++) {
			initial_cost += ComputeCost(vehicles.get(i), plans.get(i));
		}

		Candidate Some_Solution = new Candidate(vehicles, plans, taskLists, initial_cost);

// Return the generated initial candidate solution
		return Some_Solution;
	}

	private static double ComputeCost(Vehicle v, List<PD_Action> plan) {

		double cost = 0.0;

		// Follow the cities on the list of actions
		City current_city = v.getCurrentCity();

		for (PD_Action act : plan) {

			// add the cost to travel to the city
			if(act.is_pickup) {
				cost = cost + current_city.distanceTo(act.task.pickupCity) * v.costPerKm();
				current_city = act.task.pickupCity;
			}
			else {
				cost = cost + current_city.distanceTo(act.task.deliveryCity) * v.costPerKm();
				current_city = act.task.deliveryCity;
			}

		}

		return cost;
	}

	// TODO: here plan is just list of pd_actions - convert aStar plan to that?

	// TODO: calculate new Cost - should current cost be a variable? How to calculate the cost?
	private long calculatePlanCost(List<Task> tasks) {
		// TODO: should it be max?
		long cost = currentCost + distribution.weight(tasks.get(0).pickupCity, tasks.get(0).deliveryCity);
		for (Vehicle v: this.vehicles) {
			// can we do that? (TaskSet conversion)
			State currentState = generateInitialState(v, (TaskSet) tasks);
			Plan plan = aStarPlan(v, (TaskSet) tasks, currentState);
			// TODO: why null pointer and not 0
			cost = (long) plan.totalDistance()*v.costPerKm();
			// TODO: convert plan to candidate!!
			// List<List<PD_Action>> tmpPlans =
			// Candidate c = new Candidate(vehicles, plans.get(vehicles.indexOf(v)), tasks, currentCost);

		}
		// TODO: write function that computes most efficient path and therefore its cost - AStar?
		// compute best AStar solution and return its cost as "approximate" - how to include different vehicles?
		// if tasks == null return currentCost
		return cost;
	}

}
