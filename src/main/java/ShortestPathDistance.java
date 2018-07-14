import Jama.Matrix;
import org.gephi.graph.api.Edge;
import org.gephi.graph.api.Graph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by aponom on 07/12/2017.
 */
public class ShortestPathDistance implements DistanceMatrixCalculator {
    public Matrix calculateDistanceMatrix(Graph graph) {

        final HashMap <Edge, HashMap<Edge, Integer>> edgeDistanceMap = new HashMap <Edge, HashMap <Edge, Integer> >(graph.getEdgeCount());

        //calculate the distances between all pairs of edges and put its into the edgeDistanceMap

        int n = graph.getEdgeCount();
        Matrix d = new Matrix(n, n);

        for (Edge startEdge: graph.getEdges()) {

            HashMap <Edge, Integer> smallEdgeDistanceMap = new HashMap <Edge, Integer> (graph.getEdgeCount());

            List<Edge> nextStage = new ArrayList<Edge>();
            nextStage.add(startEdge);
            int stage = 0;
            smallEdgeDistanceMap.put(startEdge, stage);
            d.set(startEdge.getId() - 1, startEdge.getId() - 1, 0.0);
            while (!nextStage.isEmpty()) {
                stage++;
                List <Edge> prev = new ArrayList(nextStage);
                nextStage.clear();

                for (Edge currEdge: prev) {
                    for (Edge nextEdge: graph.getEdges(currEdge.getSource())  ) {
                        if (smallEdgeDistanceMap.get(nextEdge) == null) {
                            nextStage.add(nextEdge);
                            smallEdgeDistanceMap.put(nextEdge, stage);
                            d.set(startEdge.getId() - 1 , nextEdge.getId() - 1, stage);
                        }
                    }

                    for (Edge nextEdge: graph.getEdges(currEdge.getTarget())  ) {
                        if (smallEdgeDistanceMap.get(nextEdge) == null) {
                            nextStage.add(nextEdge);
                            smallEdgeDistanceMap.put(nextEdge, stage);
                            d.set(startEdge.getId() - 1 , nextEdge.getId() - 1 , stage);
                        }
                    }
                }
            } //while

            edgeDistanceMap.put(startEdge, smallEdgeDistanceMap);
        }

        return d;
    }

    @Override
    public String getShortName() {
        return "SP";
    }

}
