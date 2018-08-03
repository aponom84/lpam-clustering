package latna;

import Jama.Matrix;
import Jama.SingularValueDecomposition;
import org.gephi.graph.api.*;

import java.util.List;

public class Utils {

    /**
     * The difference between 1 and the smallest exactly representable number
     * greater than one. Gives an upper bound on the relative error due to
     * rounding of floating point numbers.
     */
    public static double MACHEPS = 2E-16;

    /**
     * Updates MACHEPS for the executing machine.
     */
    public static void updateMacheps() {
        MACHEPS = 1;
        do {
            MACHEPS /= 2;
        } while (1 + MACHEPS / 2 != 1);
    }

    /**
     * Computes the Moore–Penrose pseudoinverse using the SVD method.
     *
     * Modified version of the original implementation by Kim van der Linde.
     */
    public static Matrix pinv(Matrix x) {
        int rows = x.getRowDimension();
        int cols = x.getColumnDimension();
        if (rows < cols) {
            Matrix result = pinv(x.transpose());
            if (result != null) {
                result = result.transpose();
            }
            return result;
        }
        SingularValueDecomposition svdX = new SingularValueDecomposition(x);
        if (svdX.rank() < 1) {
            return null;
        }
        double[] singularValues = svdX.getSingularValues();
        double tol = Math.max(rows, cols) * singularValues[0] * MACHEPS;
        double[] singularValueReciprocals = new double[singularValues.length];
        for (int i = 0; i < singularValues.length; i++) {
            if (Math.abs(singularValues[i]) >= tol) {
                singularValueReciprocals[i] = 1.0 / singularValues[i];
            }
        }
        double[][] u = svdX.getU().getArray();
        double[][] v = svdX.getV().getArray();
        int min = Math.min(cols, u[0].length);
        double[][] inverse = new double[cols][rows];
        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < u.length; j++) {
                for (int k = 0; k < min; k++) {
                    inverse[i][j] += v[i][k] * singularValueReciprocals[k] * u[j][k];
                }
            }
        }
        return new Matrix(inverse);
    }

    public static double[][] buildLineGraph(Graph graph) {
        double a[][] = new double[graph.getEdgeCount()][graph.getEdgeCount()];

        for (Node node : graph.getNodes()) {

            for (Edge edge1 : graph.getEdges(node)) {
                for (Edge edge2 : graph.getEdges(node)) {
                    if (edge1 != edge2) {
                        a[edge1.getId() - 1][edge2.getId() - 1] = 1.0;
                        a[edge2.getId() - 1][edge1.getId() - 1] = 1.0;
                    }
                }
            }

        }

        return a;
    }

}
