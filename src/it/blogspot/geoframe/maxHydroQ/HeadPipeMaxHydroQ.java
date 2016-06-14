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

import it.blogspot.geoframe.SewerPipeDimensioning.SewerPipeDimensioning;
import it.blogspot.geoframe.hydroGeoEntities.area.DrainageArea;
import it.blogspot.geoframe.hydroGeoEntities.line.Pipe;
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

    @Description("Parameter of the recurrence interval")
    @In
    private double a;

    @Description("Parameter of the recurrence interval")
    @In
    private double n;

    private double r;
    private double n0;
    private double velocity;
    final private DrainageArea drainageArea;
    private Pipe pipe;

    /**
     * @brief Default constructor
     */
    public HeadPipeMaxHydroQ(final DrainageArea drainageArea) {
        this.drainageArea = drainageArea;
        this.pipe = this.drainageArea.getPipe();
    }

    @Execute
    public void process() {

        computeMaxFlow();

    }

    protected void initializeFirstAttemptValues() {
        pipe.setVelocity(1.0);
        r = 1.0;
    }

    protected void convergenceLoop() {
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
    }

    private double computeResidual(final double r) {
        return Math.abs(r - this.r)/this.r;
    }

    private double computeN() {
        final double product = drainageArea.getResidenceTime() * velocity;
        final double denominator = GEOunitsTransform.minutes2seconds(product);
        // @TODO: how to compute a first attempt length of the pipe? Is it
        // possibile to compute lenght just from 2D coordinates and then adjust
        // it at each loop?
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

    private double computeUdometricCoefficient(final double n0, final double r) {
        final double precipitationTime = r * drainageArea.getResidenceTime();
        return drainageArea.getUrbanRunoffCoefficient() * a * Math.pow(precipitationTime, n-2) * (1+GEOunitsTransform.minutes2seconds(CELERITYFACTOR * pipe.getVelocity() * precipitationTime/pipe.getLength()) - 1/n0 * Math.log(Math.exp(n0) + Math.exp(r) -1)) * 166.667;
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
