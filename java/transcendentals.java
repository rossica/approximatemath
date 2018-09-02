
public class transcendentals {
    /**
     * Constants for PI borrowed from
     * https://stackoverflow.com/questions/507819/pi-and-accuracy-of-a-floating-point-number
     */
    public static final double PI = 3.141592653589793D;
    public static final double PI2 = 9.869604401089357D;
    public static final double HALF_PI = 1.570796326794896D;
    public static final double THREE_HALF_PI = 4.712388980384689D;

    public static double toRadians(double deg) {
        return (deg * PI) / 180.0d;
    }

    public static double abs(double a) {
        if (a <= -0.0d) {
            return a * -1;
        }
        return a;
    }

    /***
     * Approximate sin using Bhaskara I's method according to
     * https://en.wikipedia.org/wiki/Bhaskara_I%27s_sine_approximation_formula
     * Only accurate from [0, pi]
    ***/
    public static double sin(double rad) {
        double sign = 1;
        // Massage input into the valid range.
        if (rad <= -0.0d) {
            if ((rad % (PI * 2)) > -PI) {
                sign = -1;
            }
            rad = -1 * rad;
        } else {
            if ((rad % (PI * 2)) > PI) {
                sign = -1;
            }
        }
        rad = rad % PI;

        double pi_rad = PI - rad;
        return sign * ((16 * rad * pi_rad) / (5 * PI2 - 4 * rad * pi_rad));
    }

    /***
     * Approximate cos using Bhaskara I's method according to
     * https://en.wikipedia.org/wiki/Bhaskara_I%27s_sine_approximation_formula
     * Only accurate from [-pi/2 , pi/2]
    ***/
    public static double cos(double rad) {
        // double sign = 1;
        // double rad2 = rad * rad;
        // return sign * ((PI2 - 4 * rad2) / (PI2 + rad2));
        return sin(rad + HALF_PI); // works, but feels like cheating...
    }

    /**
     * Approximate atan2 using the answer from here:
     * https://math.stackexchange.com/questions/1098487/atan2-faster-approximation
     */
    public static double atan2(double x, double y) {
        double abs_x, abs_y, a, s, r;
        if (x <= 0.0f) {
            abs_x = 0.0f - x;
        } else {
            abs_x = x;
        }

        if (y <= 0.0f) {
            abs_y = 0.0f - y;
        } else {
            abs_y = y;
        }


        if (abs_x <= abs_y) {
            a = abs_x / abs_y;
        } else {
            a = abs_y / abs_x;
        }

        s = a * a;

        r = ((-0.0464964749 * s + 0.15931422) * s - 0.327622764) * s * a + a;

        if (abs_y > abs_x) {
            // return 1.57079637 - r;
            System.out.println("first case");
            // return HALF_PI - r; // sometimes it should be this
            return  r - HALF_PI; // and sometimes it should be this.
        } else if (x < 0) {
            System.out.println("second case");
            // return PI - r; // sometimes it should be this
            return r- PI; // other times it should be this
        } else if (y < 0) {
            System.out.println("third case");
            return r * -1;
        } else {
            System.out.println("fourth case");
            return 0; // shouldn't happen. here to make javac happy.
        }
    }

    /**
     * atan2 implementation as suggested by
     * https://www.dsprelated.com/showarticle/1052.php
     * Note: this implementation actually outputs correct values, unlike above.
     */
    public static double atan(double z) {
        return (0.97239411f + -0.19194795f * z * z) * z;
    }

    public static double atan2_b(double y, double x) {
        double abs_x = (x < 0.0d) ? 0.0d - x : x;
        double abs_y = (y < 0.0d) ? 0.0d - y : y;

        if (x != 0.0d) {
            if (abs_x > abs_y) {
                final double z = y / x;
                if (x > 0.0d) {
                    // atan2(y,x) = atan(y/x) if x > 0
                    return atan(z);
                } else if (y >= 0.0d) {
                    // atan2(y,x) = atan(y/x) + PI if x < 0, y >= 0
                    return atan(z) + PI;
                } else {
                    // atan2(y,x) = atan(y/x) - PI if x < 0, y < 0
                    return atan(z) - PI;
                }
            } else { // Use property atan(y/x) = PI/2 - atan(x/y) if |y/x| > 1.
                final double z = x / y;
                if (y > 0.0d) {
                    // atan2(y,x) = PI/2 - atan(x/y) if |y/x| > 1, y > 0
                    return -atan(z) + HALF_PI;
                } else {
                    // atan2(y,x) = -PI/2 - atan(x/y) if |y/x| > 1, y < 0
                    return -atan(z) - HALF_PI;
                }
            }
        } else {
            if(y > 0.0d) { // x = 0, y > 0
                return HALF_PI;
            } else if (y < 0.0d) { // x = 0, y < 0
                return -HALF_PI;
            }
        }
        return 0.0d; // x,y = 0. Could return NaN instead.
    }

    /**
     * Fast approximation of sqrt using Newton's method as described here:
     * https://dsp.stackexchange.com/questions/17269/what-approximation-techniques-exist-for-computing-the-square-root
     * (see the answer: https://dsp.stackexchange.com/a/17314)
     */
    public static double sqrt(double a) {
        double x = 1.0d; // Initial guess. Varying this can reduce iterations.
        double x_prev; // previous iteration's guess.
        int count = 0;

        // Find a better initial guess
        while (x*x*a >= 3.0d) {
            x *= 0.5d;
        }

        do {
            x_prev = x;
            x = 0.5  * (3 * x_prev - x_prev * x_prev * x_prev * a); // Find square root of reciprocal
        } while (++count < 10); // cap based on iteration count
        // while (abs(x - x_prev) > 1e-6); // cap based desired accuracy

        return x * a; // invert the reciprocal by multiplying by a.
    }

    public static double haversine(double lon1, double lat1, double lon2, double lat2) {
        double phi1 = toRadians(lat1);
        double phi2 = toRadians(lat2);
        double deltaPhi = toRadians(lat2-lat1);
        double deltaLambda = toRadians(lon2-lon1);

        double sin2HalfDeltaPhi = sin(deltaPhi / 2);
        // System.out.println("sin(deltaPhi/2): " + sin2HalfDeltaPhi);
        sin2HalfDeltaPhi *= sin2HalfDeltaPhi;
        // System.out.println("sin^2(deltaPhi/2): " + sin2HalfDeltaPhi);

        double sin2HalfDeltaLambda = sin(deltaLambda / 2);
        // System.out.println("sin(deltaLambda/2): " + sin2HalfDeltaLambda);
        sin2HalfDeltaLambda *= sin2HalfDeltaLambda;
        // System.out.println("sin^2(deltaLambda/2): " + sin2HalfDeltaLambda);

        // System.out.println("phi1: " + phi1 + " phi2: " + phi2);

        double a = sin2HalfDeltaPhi + cos(phi1) * cos(phi2) * sin2HalfDeltaLambda;
        // System.out.println("a: " + a);
        double c = 2 * atan2_b(sqrt(a), sqrt(1-a));
        // System.out.println("c: " + c);

        return 6371 * c;
    }

    public static void main(String[] args) {
        // Test sin and cos compared to the official math library.

        // for(double i = Math.PI * -2; i <= (Math.PI * 2) + 1; i += Math.PI / 8) {
        //     System.out.println(String.format("%f:\tmath.sin = %f\tmysin = %f\t math.cos = %f\t mycos = %f",
        //     i,
        //     Math.sin(i),
        //     sin(i),
        //     Math.cos(i),
        //     cos(i)
        //     ));
        // }


        // test atan against the official math library.

        // for(double x = -1; x <= 1; ++x) {
        //     for(double y = -1; y <= 1; ++y) {
        //         System.out.println(
        //             String.format("%f, %f:\tmath.atan2 = %f\tmyatan2 = %f",
        //                 x,
        //                 y,
        //                 Math.atan2(y, x),
        //                 // atan2(x, y)
        //                 atan2_b(y, x)
        //                 ));
        //     }
        // }

        // double x = 0, y = -1;
        // System.out.println(
        //     String.format("%f, %f:\tmath.atan2 = %f\tmyatan2 = %f",
        //         x,
        //         y,
        //         Math.atan2(y, x),
        //         atan2(x, y)
        //         ));

        // x = -2;
        // y = -1;
        // System.out.println(
        //     String.format("%f, %f:\tmath.atan2 = %f\tmyatan2 = %f",
        //         x,
        //         y,
        //         Math.atan2(y, x),
        //         atan2_b(y, x)
        //         ));

        // x = 2;
        // y = -1;
        // System.out.println(
        //     String.format("%f, %f:\tmath.atan2 = %f\tmyatan2 = %f",
        //         x,
        //         y,
        //         Math.atan2(y, x),
        //         atan2(x, y)
        //         ));


        // test sqrt against the official math library

        // for(double a = 0.01d; a < 100; a += a) {
        //     System.out.println(
        //         String.format("%f:\tmath.sqrt = %f\tmysqrt = %f",
        //             a,
        //             Math.sqrt(a),
        //             sqrt(a)
        //             ));
        // }


        // test haversine

        double lat1 = 9.0;
        double lon1 = 13.0;
        double lat2 = -88.65;
        double lon2 = -122.32;
        System.out.println(
            String.format("(%f, %f) -> (%f, %f) = %fkm",
            lon1,
            lat1,
            lon2,
            lat2,
            haversine(lon1, lat1, lon2, lat2)
        ));
    }
}