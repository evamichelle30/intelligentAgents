package template;

import java.util.*;
import logist.topology.Topology.City;

public class AStar {
    public static List<City> findShortestPath(City start, City goal) {
        PriorityQueue<AStarNode> openSet = new PriorityQueue<>();
        Set<City> closedSet = new HashSet<>();
        Map<City, City> cameFrom = new HashMap<>();
        Map<City, Double> gScore = new HashMap<>();

        gScore.put(start, 0.0);
        double hScore = start.distanceTo(goal);
        openSet.add(new AStarNode(start, hScore));

        while (!openSet.isEmpty()) {
            City current = openSet.poll().city;

            if (current.equals(goal)) {
                return reconstructPath(cameFrom, current);
            }

            closedSet.add(current);

            for (City neighbor : current.neighbors()) {
                if (closedSet.contains(neighbor)) {
                    continue;
                }

                double tentativeGScore = gScore.get(current) + current.distanceTo(neighbor);

                if (!openSetContains(openSet, neighbor) || tentativeGScore < gScore.get(neighbor)) {
                    cameFrom.put(neighbor, current);
                    gScore.put(neighbor, tentativeGScore);

                    double fScore = tentativeGScore + neighbor.distanceTo(goal);
                    openSet.add(new AStarNode(neighbor, fScore));
                }
            }
        }

        return null; // No path found
    }

    private static List<City> reconstructPath(Map<City, City> cameFrom, City current) {
        List<City> path = new ArrayList<>();
        while (cameFrom.containsKey(current)) {
            path.add(current);
            current = cameFrom.get(current);
        }
        Collections.reverse(path);
        return path;
    }

    private static boolean openSetContains(PriorityQueue<AStarNode> openSet, City city) {
        return openSet.stream().anyMatch(node -> node.city.equals(city));
    }

    private static class AStarNode implements Comparable<AStarNode> {
        private final City city;
        private final double fScore;

        public AStarNode(City city, double fScore) {
            this.city = city;
            this.fScore = fScore;
        }

        @Override
        public int compareTo(AStarNode other) {
            return Double.compare(this.fScore, other.fScore);
        }
    }
}
