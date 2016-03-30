import org.jblas.DoubleMatrix;

import static org.jblas.MatrixFunctions.pow;

public class SchwarzschildBench {

    private static DoubleMatrix boltzmann(DoubleMatrix matrix) {
        return pow(matrix, 4).mul(5.67e-8).div(Math.PI);
    }


    public static void main(String[] args) {
        DoubleMatrix matrix = null;
        for(int i = 0; i <10000000; i++) {
            matrix = boltzmann(new DoubleMatrix(1, 3, 1, 2, 3));

        }
        for(Double num : matrix.toArray()) {
            System.out.println(num);
        }
    }
}
