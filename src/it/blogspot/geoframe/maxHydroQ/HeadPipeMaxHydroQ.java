/*
 * GNU GPL v3 License
 *
 * Copyright 2015 AboutHydrology (Riccardo Rigon)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package it.blogspot.geoframe.maxHydroQ;

import oms3.annotations.*;

import it.blogspot.geoframe.hydroGeoEntities.area.DrainageArea;
import it.blogspot.geoframe.hydroGeoEntities.line.Pipe;
import it.blogspot.geoframe.sewerPipeDimensioning.SewerPipeDimensioning;
import it.blogspot.geoframe.utils.GEOunitsTransform;

/**
 * @TODO: might be useful to implement a <strong>FACTORY METHOD</strong>?
 *
 * @author sidereus, francesco.serafin.3@gmail.com
 * @date May, 15th 2016
 * @copyright GNU Public License v3 GWH-2b4 (Riccardo Rigon)
 */
@Description()
// @Author()
@Keywords()
// @Label()
// @Name()
@Status()
// @License()
public class HeadPipeMaxHydroQ extends MaxHydroQ {

    private final double CELERITYFACTOR = 1;
    private final int MAXITERATION = 40;
    // @TODO: HORRIBLE HACK - Tolerance must be set as machine epsilon
    private final double TOLERANCE = 0.0001;

    @Description("Parameter of the recurrence interval")
    @In
    private double a; // [mm/min^n]

    @Description("Parameter of the recurrence interval")
    @In
    private double n; // [-]

    private double r;
    private double n0;
    final private DrainageArea drainageArea;
    private Pipe pipe;

    /**
     * @brief Default constructor
     */
    public HeadPipeMaxHydroQ(final DrainageArea drainageArea, final double a, final double n) {
        this.drainageArea = drainageArea;
        this.pipe = this.drainageArea.getPipe();
        this.a = a;
        this.n = n;
    }

    @Execute
    public void process() {

        computeMaxFlow();

    }

    protected void initializeFirstAttemptValues() {
        pipe.setVelocity(1.0);
        r = 1.0;
    }

    protected DrainageArea convergenceLoop() {
        SewerPipeDimensioning pipeDimensioning = new SewerPipeDimensioning();
        for(int iteration = 0; iteration < MAXITERATION; iteration++) {
            n0 = computeN();
            double r = computeR(n0);
            pipe.setDischarge(computeMaxDischarge(n0, r));
            pipe = pipeDimensioning.run(pipe);
            // @TODO: is it possible to check r before computing
            if (computeResidual(r) <= TOLERANCE) break;
            else this.r = r;
        }
        pipe.setPeakTime(computePeakTime());
        drainageArea.setPipe(pipe);
        return drainageArea;
    }

    /**
     * @brief Time in which the maximum between the maximum discharges is
     *        registered at the output of the designed pipe
     *
     * @return
     */
    private double computePeakTime() {
        return drainageArea.getResidenceTime() * Math.log(Math.exp(n0) + Math.exp(r) - 1);
    }

    private double computeResidual(final double r) {
        return Math.abs(r - this.r)/this.r;
    }

    private double computeN() {
        final double residenceTime = GEOunitsTransform.minutes2seconds(drainageArea.getResidenceTime());
        final double denominator = residenceTime * pipe.getVelocity();
        return pipe.getLength() / denominator;
    }

    private double computeR(final double n0) {
        final double n0_inf = 0.1;
        final double n0_sup = 2 * n0 + 5;
        return bisection(n0, n0_inf, n0_sup);
    }

    private double computeMaxDischarge(final double n0, final double r) {
        return computeUdometricCoefficient(n0, r) * drainageArea.getArea();
    }

    /**
     * @TODO: might the computing of the udometric coefficient be moved in the
     * drainage area class?
     *
     * @param n0
     * @param r
     * @return
     */
    private double computeUdometricCoefficient(final double n0, final double r) {
        final double precipitationTime = r * drainageArea.getResidenceTime();
        final double udometricCoefficient = drainageArea.getUrbanRunoffCoefficient() * a * Math.pow(precipitationTime, n-1) * (1+CELERITYFACTOR * pipe.getVelocity() * GEOunitsTransform.minutes2seconds(precipitationTime)/pipe.getLength() - 1/n0 * Math.log(Math.exp(n0) + Math.exp(r) -1));
        return GEOunitsTransform.cubicMeters2liters(GEOunitsTransform.hectars2meters(GEOunitsTransform.millimiters2meters(GEOunitsTransform.seconds2minutes(udometricCoefficient))));
    }

    /**
     * @TODO: put bisection method inside GEOframeUtils package. To implement it,
     * using the java 8 feauture to pass functions
     *
     * @param n0
     * @param n0_inf
     * @param n0_sup
     * @return
     */
    private double bisection(final double n0, final double n0_inf, final double n0_sup) {
            final double function = distributedInflux(n0, n0_inf);
            double function_mid = distributedInflux(n0, n0_sup);

            double dx;
            double rtb = (function < 0) ? n0_inf : n0_sup;
            dx = (function < 0) ? (n0_sup - n0_inf) : (n0_inf - n0_sup);

            final int JMAX = 40;
        try {
            for (int j=0; j < JMAX; j++) {
                dx *= 0.5;
                final double xmid = rtb + dx;
                function_mid = distributedInflux(n0, xmid);

                if (function_mid <= 0) rtb = xmid;
                if (Math.abs(dx) < TOLERANCE || function_mid == 0) break;
            }
        } catch(RuntimeException exception) {
            throw new RuntimeException(exception.getMessage());
        }

        return rtb;
    }

    private double distributedInflux(final double n0, final double r) {
        final double numerator = r * (1 - Math.exp(r) / (Math.exp(n0) + Math.exp(r) - 1));
        final double denominator = n0 + r - Math.log(Math.exp(n0) + Math.exp(r) - 1);

        return 1-n-numerator / denominator;
    }

}
