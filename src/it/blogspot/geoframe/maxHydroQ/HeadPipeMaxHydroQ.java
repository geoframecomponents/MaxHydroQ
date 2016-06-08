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

import it.blogspot.geoframe.utils.GEOunitsTransform;

/**
 * @mainpage Maximum flow from Hydrograph
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
    private final double TOLERANCE = 0.0001;

    @Description("Mean residence time")
    @In
    private double residenceTime;

    @Description("Length of the pipe")
    @In
    private double pipeLength;

    @Description("Velocity of the water in the drainage area")
    @In
    private double udometricCoefficient;

    @Description("Influx coefficient to the net")
    @In
    private double influxCoefficient;

    @Description("Parameter of the recurrence interval")
    @In
    private double a;

    @Description("Parameter of the recurrence interval")
    @In
    private double n;

    @Description("Drainage area")
    @In
    private double drainageArea;

    @Description("Maximum output discharge from the drainage area")
    @Out
    private double maxDischarge;

    private double r;
    private double n0;
    private double diameter;
    private double velocity;
    /**
     * @brief Default constructor
     */
    HeadPipeMaxHydroQ() {}

    @Execute
    public void process() {

        computeMaxFlow();

    }

    protected void initializeFirstAttemptValues() {
        diameter = 0.0;
        velocity = 1.0;
        r = 1.0;
    }

    protected void convergenceLoop() {
        for(int iteration = 0; iteration < MAXITERATION; iteration++) {
            n0 = computeN();
            double r = computeR(n0);
            maxDischarge = computeMaxDischarge(n0, r);
            // compute diameter
            if (computeResidual(r) <= TOLERANCE) break;
            else this.r = r;
        }
    }

    private double computeResidual(final double r) {
        return Math.abs(r - this.r)/this.r;
    }

    private double computeN() {
        return pipeLength / GEOunitsTransform.minutes2seconds(residenceTime * velocity);
    }

    private double computeR(final double n0) {
        final double n0_inf = 0.1;
        final double n0_sup = 2 * n0 + 5;
        return bisection(n0, n0_inf, n0_sup);
    }

    private double computeMaxDischarge(final double n0, final double r) {
        return computeUdometricCoefficient(n0, r) * drainageArea;
    }

    private double computeUdometricCoefficient(final double n0, final double r) {
        final double precipitationTime = r * residenceTime;
        return influxCoefficient * a * Math.pow(precipitationTime, n-1) * (1+GEOunitsTransform.minutes2seconds(CELERITYFACTOR * velocity * precipitationTime/pipeLength) - 1/n0 * Math.log(Math.exp(n0) + Math.exp(r) -1)) * 166.667;
    }

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
        return 1-n-numerator(n0, r)/denominator(n0, r);
    }

    private double numerator(final double n0, final double r) {
        return r * (1 - Math.exp(r) / (Math.exp(n0) + Math.exp(r) - 1));
    }

    private double denominator(final double n0, final double r) {
        return n0 + r - Math.log(Math.exp(n0) + Math.exp(r) - 1);
    }

}
