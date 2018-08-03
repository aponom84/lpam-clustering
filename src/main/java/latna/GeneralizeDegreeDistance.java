package latna;

import Jama.Matrix;
import org.gephi.graph.api.Edge;
import org.gephi.graph.api.Graph;

import java.util.HashMap;

/**
 * Created by aponom on 07/12/2017.
 */
public class GeneralizeDegreeDistance implements DistanceMatrixCalculator {
    @Override
    public Matrix calculateDistanceMatrix(Graph graph) {
            final int maxDiameter = 20;
            System.out.println("calculating GeneralizeDegreeDistanceMatrix");

            //calculate line graph
            final HashMap <Edge, Integer> edgeMap = new HashMap();

            //fill the edge map
            for (Edge edge: graph.getEdges())
                edgeMap.put(edge, edge.getId());


            //--------------------start of distance matrix calculation ---------------------------------
            double a[][] = Utils.buildLineGraph(graph);

            int n = a[0].length;

            Matrix A = new Matrix(a);
            Matrix AA = new Matrix(a);

            long ak [][][] = new long [n][n][maxDiameter];
            double gd[][] = new double [n][maxDiameter];
            double mk[] = new double [maxDiameter];
            long sp[][] = new long [n][n];

            int diameter  = 1;

            boolean needToContinue;

            do {
                //A.print(writer, 3, 1);  writer.flush();
                needToContinue = false;

                for (int i = 0; i < n; i++)
                    for (int j = i+1; j < n; j++ ) {
                        if (i == j) continue;

                        if ( (sp[i][j] == 0) &&  (A.get(i, j) != 0.0) ) {
                            needToContinue = true;
                            sp[i][j] = diameter;
                            sp[j][i] = diameter;
                            ak[i][j][diameter] = (long) A.get(i, j);
                            ak[j][i][diameter] = (long) A.get(j, i);

                            gd[i][diameter]+= ak[i][j][diameter];
                            gd[j][diameter]+= ak[j][i][diameter];
                            mk[diameter] += 2*ak[i][j][diameter];
                        }
                    }

                A = A.times(AA);
                diameter++;
            } while (needToContinue);

            diameter--;
            diameter--;

            System.out.println("Diameter: " + diameter);

            Matrix D = new Matrix(n, n);
            for (int i = 0; i < n; i++) {
                D.set(i, i, 0.0);
                for (int j = i+1; j < n; j++  ) {
                    double distance  = 0;
/*
                    for (int k = 1; k <= diameter; k++ ) {
                        distance +=   k * ( (double) gd[i][k] ) * ( (double) gd[j][k] / mk[k] );
                    }
*/

                   for (int k = 1; k <= diameter; k++ ) {
                       distance += k * ( (double) gd[i][k] / mk[k] ) * ( (double) gd[j][k] / mk[k] );
                   }
/*
                    for (int k = 1; k <= sp[i][j]; k++ ) {
                        distance += k * ( (double) gd[i][k] / mk[k] ) * ( (double) gd[j][k] / mk[k] );
                    }
                    */
//                      for (int k = (int) sp[i][j]; k <= diameter; k++ ) {
//                           distance += k * ( (double) gd[i][k] / mk[k] ) * ( (double) gd[j][k] / mk[k] );
//                      }

                    if (distance == 0.0) {
                        System.out.println("Error!!!");
                        System.exit(2);
                    }
                    D.set(i, j, distance); D.set(j, i, distance);
                }
            }

            System.out.println("edgeDistanceMap has been calculated\n");

            return D;

        }

    @Override
    public String getShortName() {
        return "GD";
    }
}
