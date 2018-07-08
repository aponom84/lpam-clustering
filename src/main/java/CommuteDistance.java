import Jama.Matrix;
import org.gephi.graph.api.Edge;
import org.gephi.graph.api.Graph;

import java.util.HashMap;

/**
 * Created by aponom on 07/12/2017.
 */
public class CommuteDistance implements DistanceMatrixCalculator {
    public Matrix calculateDistanceMatrix(Graph graph) {
        final HashMap <Edge, HashMap<Edge, Double>> edgeDistanceMap = new HashMap <Edge, HashMap <Edge, Double> >(graph.getEdgeCount());

        //claculate lineral graph
        HashMap <Edge, Integer> edgeMap = new HashMap();

        //fill the edge map
        for (Edge edge: graph.getEdges())
            edgeMap.put(edge, edge.getId());

        //--------------------start of distance maxtrix calculation ---------------------------------
        double a[][] = Utils.buildLineGraph(graph);

        int n = graph.getEdgeCount(); //now we have

        double V = 0;
        Matrix D = new Matrix(n, n);
        for (int i=0; i < n; i++) {
            edgeDistanceMap.put(graph.getEdge(i+1), new HashMap());
            double degree = 0;
            for (int j=0; j < n; j++  ) {
                degree = degree+a[i][j];
                V = V + a[i][j];
            }
            D.set(i, i, degree);
            System.out.println(degree);
        }

        Matrix A = new Matrix(a);
        Matrix L = D.minus(A);

        L = Utils.pinv(L);

        Matrix d = new Matrix(n, n);

        for (int i = 0; i < n; i++) {

            Matrix ei = new Matrix(1,n);
            ei.set(0, i, 1);
            for (int j = 0; j < n; j++  ) {
                Matrix ej = new Matrix(1,n);
                ej.set(0, j, 1);

                Matrix left = ei.minus(ej);

                left = left.times(V);

                left =  left.times(L);

                Matrix right = ei.transpose().minus(ej.transpose());

                left = left.times(right);

                d.set(i, j, Math.sqrt(left.get(0, 0)));
                edgeDistanceMap.get(graph.getEdge(i+1)).put(graph.getEdge(j+1),    Math.sqrt(left.get(0, 0))  );
                edgeDistanceMap.get(graph.getEdge(j+1)).put(graph.getEdge(i+1),    Math.sqrt(left.get(0, 0))  );
            }
        }
        return d;
    }

    @Override
    public String getShortName() {
        return "CM";
    }
}
