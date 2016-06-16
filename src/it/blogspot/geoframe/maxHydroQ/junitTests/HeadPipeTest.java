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
package it.blogspot.geoframe.maxHydroQ.junitTests;

import org.junit.Test;

import it.blogspot.geoframe.hydroGeoEntities.area.DrainageArea;
import it.blogspot.geoframe.hydroGeoEntities.line.Pipe;
import it.blogspot.geoframe.hydroGeoEntities.point.HydroGeoPoint;
import it.blogspot.geoframe.hydroGeoEntities.point.InspectionChamber;
import it.blogspot.geoframe.maxHydroQ.HeadPipeMaxHydroQ;
import it.blogspot.geoframe.maxHydroQ.MaxHydroQ;

/**
 *
 *
 * @author sidereus, francesco.serafin.3@gmail.com
 * @version 0.1
 * @date June 16, 2016
 * @copyright GNU Public License v3 GWH-2b4
 */
public class HeadPipeTest {

    @Test
    public void computeHeadPipe() {
        final double ks = 65;
        final double fillCoefficient = 0.80;
        final HydroGeoPoint startInspectionChamber = new InspectionChamber(10,10,0);
        final HydroGeoPoint endInspectionChamber = new InspectionChamber(130,10,0);

        Pipe pipe = new Pipe(ks, fillCoefficient, startInspectionChamber, endInspectionChamber);
        final double area = 1.937;
        final double urbanRunoffCoefficient = 0.7;
        final double alpha = 0.39;
        final double averageSlope = 1.2;
        final double a = 4.97;
        final double n = 0.61;
        DrainageArea drainageArea = new DrainageArea(pipe, area, urbanRunoffCoefficient, alpha, averageSlope);

        MaxHydroQ pipeDesign = new HeadPipeMaxHydroQ(drainageArea, a, n);
        drainageArea = pipeDesign.computeMaxFlow();

        System.out.println("Pipe diameter: " + drainageArea.getPipe().getDiameter());
        System.out.println("Pipe slope: " + drainageArea.getPipe().getSlope());
        System.out.println("Pipe velocity: " + drainageArea.getPipe().getVelocity());
        System.out.println("Residence time: " + drainageArea.getResidenceTime());
    }

}
