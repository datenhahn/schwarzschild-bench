import org.apache.commons.lang3.time.StopWatch;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;

import static org.jblas.MatrixFunctions.pow;

public class SchwarzschildBench {

    private static final double GRAVITY = 9.81;
    private static final double CP = 1.006e3;


//    #@vectorize([float64(float64, float64, float64)])
//    def T2theta(T, p, p0):
//            """ convert Temperature in K to potential Temperature at reference pressure p0 """
//            return T * (p0/p)**(2./7)

    private static DoubleMatrix bigT2theta(DoubleMatrix bigT, DoubleMatrix p, double p0) {
        return bigT.mul(MatrixFunctions.pow((DoubleMatrix.scalar(p0).div(p)), (2. / 7)));
    }

    private static DoubleMatrix theta2BigT(DoubleMatrix bigT, DoubleMatrix p, double p0) {
        return bigT.div(MatrixFunctions.pow((DoubleMatrix.scalar(p0).div(p)), (2. / 7)));
    }

    private static double rmse(DoubleMatrix a, DoubleMatrix b) {
        return Math.sqrt(MatrixFunctions.pow(a.sub(b), 2).mean());
    }
//            #@jit('float64(float64[:], float64[:])', nopython=True)
//    def rmse(a, b):
//            return np.sqrt(np.mean((a-b)**2))

    private static DoubleMatrix stephanboltzmann(DoubleMatrix matrix) {
        return pow(matrix, 4).mul(5.67e-8).div(Math.PI);
    }

    /**
     * def schwarzschild(dtau, w0, g, B, Ag, Eup, Edn, divE):
     *
     * Compute radiative transfer according to schwarzschild equation, neglecting scattering
     * returns Irradiace downward, upward and divergence of radiation
     *
     * @return irradianceDivergency
     */
    private static void schwarzschild(DoubleMatrix dTau, DoubleMatrix w0, DoubleMatrix g, DoubleMatrix bigB,
                                      int ag, DoubleMatrix eUp, DoubleMatrix eDown, DoubleMatrix divergencyE)
    {

        double nmu = 5;
        double dmu = 1 / nmu;

        int ke = dTau.length;

        DoubleMatrix kabs = dTau.mul(DoubleMatrix.scalar(1).sub(w0));

        for (int imu = 0; imu < nmu; imu++) {
            double mu = (imu + 0.5) * dmu;
            DoubleMatrix bigT = MatrixFunctions.exp(kabs.mul(-1).div(mu));
            double lUp = bigB.get(ke - 1) * (1 - ag);
            eUp.put(ke, eUp.get(ke) + lUp * mu);

            double lDn = 0;
            eDown.put(0, eDown.get(0) + lDn * mu);

            for (int k = ke - 1; k > -1; k--) {
                lUp = lUp * bigT.get(k) + bigB.get(k) * (1 - bigT.get(k));
                eUp.put(k, eUp.get(k) + lUp * mu);
            }

            for (int k = 0; k < ke; k++) {
                lDn = lDn * bigT.get(k) + bigB.get(k) * (1 - bigT.get(k));
                eDown.put(k + 1, eDown.get(k + 1) + lDn * mu);
            }
        }

        eUp.muli(2).muli(Math.PI).muli(dmu);
        eDown.muli(2).muli(Math.PI).muli(dmu);
        divergencyE.copy(
                eUp.getRange(1, eUp.rows)
                   .sub(eUp.getRange(0, eUp.rows - 1))
                   .add(eDown.getRange(0, eDown.rows - 1)).sub(eDown.getRange(1, eDown.rows)));
        divergencyE.put(divergencyE.rows-1, divergencyE.get(divergencyE.rows - 1) + eDown.get(eDown.rows                                                                                              - 1)
                                            - eUp.get(eUp.rows - 1));
    }
//        for imu in xrange(Nmu):
//            mu = (imu + .5)*dmu
//
//            T = np.exp(-kabs / mu)
//
//            Lup = B[ke-1]*(1.-Ag)
//            Eup[ke] = Eup[ke] + Lup*mu
//
//            Ldn = 0
//                Edn[0] = Edn[0] + Ldn*mu
//
//            for k in xrange(ke-1, -1, -1):
//                Lup = Lup * T[k] + B[k] * (1.-T[k])
//                Eup[k] = Eup[k] + Lup*mu
//
//            for k in xrange(0, ke):
//                Ldn = Ldn * T[k] + B[k] * (1.-T[k])
//                Edn[k+1] = Edn[k+1] + Ldn*mu
//
//        Eup[:] = Eup[:] * 2 * np.pi * dmu
//        Edn[:] = Edn[:] * 2 * np.pi * dmu
//        divE[:] = Eup[1:] - Eup[:-1] + Edn[:-1] - Edn[1:]
//        divE[-1] += Edn[-1] - Eup[-1]

    private static void cmfj() {
        cmfj(100, 288, 101300);
    }

    private static void cmfj(int nlyr, int t0, int p0) {

        StopWatch watch = new StopWatch();
        watch.start();
        DoubleMatrix pt = DoubleMatrix.linspace(0, p0, nlyr + 1);
        DoubleMatrix pm = pt.getRange(1, pt.rows).add(pt.getRange(0, pt.rows - 1)).div(2);

        DoubleMatrix bigT = new DoubleMatrix(nlyr);

        int index = 0;
        for (int i = nlyr - 1; i >= 0; i--) {
            int value = i * -1;
            bigT.put(index, value);
            index++;
        }
        bigT.addi(t0);
        DoubleMatrix lastBigT = new DoubleMatrix(nlyr);
        lastBigT.copy(bigT);

        DoubleMatrix divE = DoubleMatrix.zeros(nlyr);

        boolean converged = false;

        while (!converged) {
            DoubleMatrix dTau = DoubleMatrix.ones(nlyr).mul(1).div(nlyr);
            DoubleMatrix w0 = DoubleMatrix.zeros(nlyr);
            DoubleMatrix g = DoubleMatrix.zeros(nlyr);
            DoubleMatrix b = stephanboltzmann(bigT);


            DoubleMatrix eUp = DoubleMatrix.zeros(nlyr + 1);
            DoubleMatrix eDown = DoubleMatrix.zeros(nlyr + 1);

            schwarzschild(dTau, w0, g, b, 0, eUp, eDown, divE);

            divE.put(divE.rows-1, divE.get(divE.rows-1) + 235);

            DoubleMatrix
                    tInc =
                    divE.mul(GRAVITY).div(CP).div(pt.getRange(1, pt.rows).sub(pt.getRange(0, pt.rows - 1)));
            double dt = 1. / MatrixFunctions.abs(tInc).max();

            bigT.addi(tInc.mul(dt));

            DoubleMatrix new_theta = reverse(bigT2theta(bigT, pm, p0).sort());
            bigT = theta2BigT(new_theta, pm, p0);

            double residual = rmse(lastBigT, bigT);
            System.out.println("residual: " + residual);

            if (residual < 1e-4) {
                converged = true;
            }
            lastBigT.copy(bigT);
        }
        watch.stop();
        System.out.println("== Converged ==");

        for (int k = 0; k < nlyr; k++) {
            System.out.println("Pressure " + pm.get(k) + "  Temperature " + bigT.get(k));
        }

        System.out.println("Time: " + watch.getTime() + "ms");
        //        print "== Converged =="
//        for k in xrange(nlyr):
//        print k, " Pressure %s  Temperature %s" % (pm[k], T[k])
//
//        print "== Time: "+ str(delta.total_seconds() *1000)+ " ms"


    }

    public static DoubleMatrix reverse(DoubleMatrix input) {
        if (input.columns > 1) {
            throw new IllegalArgumentException("only supports 1 column matrices");
        } else {
            DoubleMatrix result = new DoubleMatrix(input.rows);
            for (int i = 0; i < input.rows; i++) {
                result.put(i, input.get(input.rows - i - 1));
            }
            return result;
        }
    }

//        while not converged:
    //        dtau = 1. * np.ones(nlyr) / nlyr
    //        w0 = np.zeros(nlyr)
    //        g = np.zeros(nlyr)
    //        B = stephanboltzmann(T)
    //
    //        Eup = np.zeros(nlyr+1)
    //        Edn = np.zeros(nlyr+1)
    //        divE = np.zeros(nlyr)
    //
    //        schwarzschild(dtau, w0, g, B, 0, Eup, Edn, divE)
    //        #        jitschwarzschild(dtau, w0, g, B, np.float64(0), Eup, Edn, divE)
    //
    //        divE[-1] += 235
    //
    //        Tinc = divE *CONST_G / CONST_CP / (pt[1:]-pt[:-1])
    //
    //        dt = 1. / np.max(np.abs(Tinc))
    //
    //        T = T + Tinc*dt
    //
    //        new_theta = np.sort(T2theta(T, pm, p0))[::-1]
    //        T = theta2T(new_theta, pm, p0)
    //
    //        residual = rmse(lastT, T)
    //        print residual
    //        if residual < 1e-4:
    //        converged = True
    //        lastT = T.copy()
//
//        b = datetime.datetime.now()
//
//        delta = b - a
//
//        print "== Converged =="
//        for k in xrange(nlyr):
//        print k, " Pressure %s  Temperature %s" % (pm[k], T[k])
//
//        print "== Time: "+ str(delta.total_seconds() *1000)+ " ms"


    public static void main(String[] args) {
//

//        DoubleMatrix divergencyE = DoubleMatrix.zeros(numberOfLayers);
//        DoubleMatrix matrix = null;
//        for (int i = 0; i < 10000000; i++) {
//            matrix = stephanboltzmann(new DoubleMatrix(1, 3, 1, 2, 3));
//
//        }
//        for (Double num : matrix.toArray()) {
//            System.out.println(num);
//        }
        cmfj();
    }
}
