package biz.personalAcademics.ellipsoid;

import java.util.ArrayList;
import biz.personalAcademics.ellipsoid.customExceptions.*;

public class EllipsoidTinyVolumes extends EllipsoidEvenDistributionRectCoords{
	
	public EllipsoidTinyVolumes(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisStart, double radianMeasureOffZAxisEnd,
			double a, double b, double c) {
		
		super(startRadianTheta, endRadianTheta, radianMeasureOffZAxisStart, radianMeasureOffZAxisEnd,
			a, b, c);
	}
	
	/**
	 * Uses a uniformly distributed set of 3 dimensional coordinates to estimate
	 * volume. Formula for number of points is = [(number of sample points on 1
	 * positive axis)^3] * 8. Due to computation time, samples are limited to
	 * ~2,700,000 sample points.
	 * 
	 * Method still under construction
	 * 
	 * @return
	 */
	public double getEstimatedVolume() {

		ArrayList<double[]> coord = generateUniformDistributionOfSamplePoints();
		
		double sampleSize = coord.size();

		double rectangularPrisimVolume = this.a / sampleSize * this.b
				/ sampleSize * this.c / sampleSize;
		double shapeVolume = 0;

		for (int i = 0; i < coord.size(); i++) {
			try {
				if (determineIfPointIsInsideShapeNonRandomDistribution(coord
						.get(i))) {
					shapeVolume += rectangularPrisimVolume;
				}
			} catch (PointOnEdgeException | PositiveZAxisException e) {
				shapeVolume += rectangularPrisimVolume / 2;
			}

		}

		// portion of volume of 3D system
		return shapeVolume;
	}

}
