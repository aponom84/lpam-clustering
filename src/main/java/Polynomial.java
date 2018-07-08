

import java.util.HashMap;
import java.util.SortedSet;

public class Polynomial {
    HashMap<SortedSet<Integer>,Double> coefficients;
    Polynomial()
    {
        coefficients = new HashMap<SortedSet<Integer>,Double>();
    }

    HashMap<SortedSet<Integer>,Double> getCoefficients()
    {
        return coefficients;
    }
    void putCoef(SortedSet<Integer> var, Double coef) {
        if (coef != 0.0) coefficients.put(var,coef);
    }
    double getCoef(SortedSet<Integer> var) {
        return coefficients.get(var);
    }
    void plus(Polynomial other) {
        for ( SortedSet<Integer> key : other.coefficients.keySet() ) {
            double val = other.coefficients.get(key);
            if( this.coefficients.get(key) == null)
                this.putCoef(key, val);
            else this.putCoef(key, this.getCoef(key) + val);
        }
    }
    void print() {
        int count = 0;
        int size = coefficients.keySet().size();
        for ( SortedSet<Integer> key : coefficients.keySet() )
        {
            count++;
            String beautyPoly = "";
            for (Integer y: key) {
                if (y != 0) {
                    beautyPoly += "y" + String.valueOf(y);
                }
            }
            System.out.print(coefficients.get(key) + beautyPoly + (count==size ? "":" + "));
        }
        System.out.println();
    }
}
