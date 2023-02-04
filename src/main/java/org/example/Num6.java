package org.example;

import jeigen.DenseMatrix;
import smile.math.matrix.Matrix;

import static jeigen.Shortcuts.zeros;

public class Num6 {

    public static void main(String[] args){
        System.out.println();
        System.out.println("a)");
        System.out.println();
        Check.calcEigenValues(getMatrixM());
        System.out.println();
        Calc.calcEigenValues(getMatrixM());
        System.out.println();
        System.out.println();
        System.out.println("b)");
        System.out.println();
        Check.calcEigenValuesAndVectors(getMatrixB());
        System.out.println();
        Calc.calcEigenValuesAndVectors(getMatrixB());
    }

    private static double[][] getMatrixM(){
        return new double[][] {{3, 6, 6, 9},
                               {1, 4, 0, 9},
                               {0, 0.2, 6, 12},
                               {0, 0, 0.1, 6}};
    }

    private static double[][] getMatrixB(){
        return new double[][] {{3, 4, 2, 4},
                               {4, 7, 1, -3},
                               {2, 1, 3, 2},
                               {4, -3, 2, 2}};
    }
}

class Calc {

    private static final double EPSILON = 1e-10;
    
    public static void calcEigenValues(double[][] matrix){
        try {

            int i = 1000;
            while(!isDiagonalMatrix(matrix)){
                Matrix smileMatrix = new Matrix(matrix);
                Matrix.QR smileQR = smileMatrix.qr();
                double[][] Q = smileQR.Q().toArray();
                double[][] R = smileQR.R().toArray();
                matrix = multiplyMatrices(R, Q);
                if(i == 0)
                    break;
                i--;
            }
        } catch (java.lang.UnsatisfiedLinkError e){
            System.out.println(e);
            return;
        }

        System.out.println();
        System.out.println("Wartosci wlasne obliczone:");
        System.out.println(matrix[0][0] + ", " + matrix[1][1]
                            + ", " + matrix[2][2] + ", " + matrix[3][3]);
    }

    public static void calcEigenValuesAndVectors(double[][] matrix){
        double[] b = new double[] {1, 1,1 ,1};
        double valueBiggest = 0;
        double valueBiggest2 = 1000;
        double eigenVal = 0;
        while(valueBiggest - valueBiggest2 > EPSILON || valueBiggest - valueBiggest2 < EPSILON *-1){
            valueBiggest = valueBiggest2;
            double[] ab = multiply(matrix, b);
            double normAb = calculateNorm(ab);
            eigenVal = normAb;
            b = divide(ab, normAb);
            valueBiggest2 = dotProduct(b, new double[] {1,1,1,1});
        }
        System.out.println();
        System.out.println("najwieksza wartosc wlasna:");
        System.out.println(eigenVal);
        System.out.println();
        System.out.println("Wektor wlasny najwiekszej wartosci:");
        PrintStuff.printVerticalVector(b);
    }

    private static boolean isDiagonalMatrix(double[][] matrix){
        boolean isUpperDiagonal = true;
        for(int i = 0; i < matrix.length-1; i++){
            if(matrix[i][i+1] > EPSILON || matrix[i][i+1] < -1 * EPSILON){
                isUpperDiagonal = false;
                break;
            }
        }
        if(isUpperDiagonal)
            return true;

        for(int i = 0; i < matrix.length-1; i++){
            if(matrix[i][i+1] > EPSILON || matrix[i][i+1] < -1 * EPSILON)
                return false;
        }
        return true;
    }

    
    private static double[][] multiplyMatrices(double[][] matrixA, double[][] matrixB) {
        int aRows = matrixA.length;
        int aColumns = matrixA[0].length;
        int bRows = matrixB.length;
        int bColumns = matrixB[0].length;

        double[][] result = new double[aRows][bColumns];

        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bColumns; j++) {
                double sum = 0;
                for (int k = 0; k < aColumns; k++) {
                    sum += matrixA[i][k] * matrixB[k][j];
                }
                result[i][j] = sum;
            }
        }
        return result;
    }

    private static double[] multiply(double[][] matrix, double[] vector) {
      int rows = matrix.length;
      int columns = matrix[0].length;
      double[] result = new double[rows];

      for (int i = 0; i < rows; i++) {
         double sum = 0;
         for (int j = 0; j < columns; j++) {
            sum += matrix[i][j] * vector[j];
         }
         result[i] = sum;
      }

      return result;
   }

   private static double calculateNorm(double[] vector) {
      double sum = 0;
      for (double v : vector) {
         sum += v * v;
      }
      return Math.sqrt(sum);
   }

   private static double[] divide(double[] vector, double scalar) {
       double[] d = new double [vector.length];
      for (int i = 0; i < vector.length; i++) {
         d[i] = vector[i] / scalar;
      }
      return d;
   }

   public static double dotProduct(double[] vector1, double[] vector2) {
    int length = Math.min(vector1.length, vector2.length);
    double result = 0;
    for (int i = 0; i < length; i++) {
        result += vector1[i] * vector2[i];
    }
    return result;
}
}

class Check {

    public static DenseMatrix.EigenResult calcEigenValues(double[][] matrix){
        DenseMatrix dMatrix = Convertor.getDenseMatrix(matrix);
        DenseMatrix.EigenResult result = dMatrix.eig();
        System.out.println();
        System.out.println("Wartosci wlasne obliczone przez biblioteke Jeigen:");
        System.out.println(result.values.getReal(0,0) + ", " + result.values.getReal(1, 0)
                            + ", " + result.values.getReal(2, 0) + ", " + result.values.getReal(3, 0));
        return result;
    }
    
    public static void calcEigenValuesAndVectors(double[][] matrix){
        DenseMatrix.EigenResult result = calcEigenValues(matrix);
        System.out.println();
        System.out.println("najwieksza wartosc wlasna obliczone przez biblioteke Jeigen:");
        int biggest = biggest(result);
        System.out.println(result.values.getReal(biggest,0));
        System.out.println();
        System.out.println("Wektor wlasny najwiekszej wartosci obliczone przez biblioteke Jeigen:");
        double[] vectors = Convertor.getMatrix(result.vectors.real(), biggest);
        PrintStuff.printVerticalVector(vectors);

    }

    private static int biggest(DenseMatrix.EigenResult res){
        int big = 0;
        double val = plusValue(res.values.getReal(0,0));
        for(int i = 1; i < 4; i++){
            double nextVal = plusValue(res.values.getReal(i,0));
            if(val < nextVal){
                val = nextVal;
                big = i;
            }
        }
        return big;
    }

    private static double plusValue(double v){
        if (v < 0 )
            return v * -1;
        return v;
    }
}

class Convertor {

    public static DenseMatrix getDenseMatrix(double[][] doubles){
        String convertor = "";
        for(double[] row : doubles){
            for(double d : row)
                convertor += d + " ";
            convertor.trim();
            convertor += "; ";
        }
        convertor = convertor.substring(0, convertor.length()-3);

        return new DenseMatrix(convertor);
    }

    public static DenseMatrix getDenseMatrixFromVerticalVector(double[] doubles){
        String convertor = "";
        for(double d : doubles)
            convertor += d + "; ";
        convertor = convertor.substring(0, convertor.length()-3);
        return new DenseMatrix(convertor);
    }

    public static DenseMatrix getDenseMatrixFromStripMatrix(double[][] doubles){
        DenseMatrix matrix = zeros(doubles[0].length, doubles[0].length);
        for(int i = 0; i < doubles[0].length; i++){
            if(i < doubles[0].length-1)
                matrix.set(0+i, 1+i, doubles[0][i]);
            matrix.set(0+i, 0+i, doubles[1][i]);
        }
        return matrix;
    }

    public static double[] getMatrix(DenseMatrix dMatrix, int i){
        double[] d = new double[4];
            for(int j = 0; j < 4; j++){
                d[j] = dMatrix.get(j, i);
            }
        return d;
    }

    public static double[] getVector(DenseMatrix dMatrix){
        return dMatrix.getValues();
    }

}

class PrintStuff {

    public static void printMatrix(double[][] doubles){
        for(double[] row : doubles) {
            printHorizontalVector(row);
            System.out.println();
        }
    }

    public static void printMatrix(float[][] doubles){
        for(float[] row : doubles) {
            printHorizontalVector(row);
            System.out.println();
        }
    }

    public static void printMatrix(int[][] doubles){
        for(int[] row : doubles) {
            printHorizontalVector(row);
            System.out.println();
        }
    }

    public static void printVerticalVector(double[] doubles){
        for(double d : doubles)
            System.out.println(d);
    }

    public static void printVerticalVector(float[] doubles){
        for(float d : doubles)
            System.out.println(d);
    }

    public static void printVerticalVector(int[] doubles){
        for(int d : doubles)
            System.out.println(d);
    }

    public static void printHorizontalVector(double[] doubles){
        for(double d : doubles)
            System.out.print(d + "\t");
    }

    public static void printHorizontalVector(float[] doubles){
        for(float d : doubles)
            System.out.print(d + "\t");
    }

    public static void printHorizontalVector(int[] doubles){
        for(int d : doubles)
            System.out.print(d + "\t");
    }
}