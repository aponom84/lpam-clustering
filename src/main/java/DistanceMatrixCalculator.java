import Jama.Matrix;
import org.gephi.graph.api.Graph;

/**
 * Created by aponom on 07/12/2017.
 */
public interface DistanceMatrixCalculator {
    public Matrix calculateDistanceMatrix(Graph graph);
    public String getShortName();
}
