package biz.personalAcademics.ellipsoid;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import biz.personalAcademics.ellipsoid.customExceptions.*;

public class EllipsoidEvenDistributionRectCoords extends EllipsoidRectangularCoords {
	
	public EllipsoidEvenDistributionRectCoords(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisStart, double radianMeasureOffZAxisEnd,
			double a, double b, double c) {
		
		super(startRadianTheta, endRadianTheta, radianMeasureOffZAxisStart, radianMeasureOffZAxisEnd,
			a, b, c);
	}
	
	/**
	 * Uses a uniformly distributed set of 3 dimentional coordinates to estimate
	 * volume. Formula for number of points is = [(number of sample points on 1
	 * positive axis)^3] * 8. Due to computation time, samples are limited to
	 * ~2,700,000 sample points.
	 * 
	 * @return
	 */
	public double getEstimatedVolume() {
		double insideShape = 0;
		ArrayList<double[]> coord = generateUniformDistributionOfSamplePoints();

		for (int i = 0; i < coord.size(); i++) {
			try {

				if (determineIfPointIsInsideShapeNonRandomDistribution(coord
						.get(i))) {
					insideShape += 1;
				}

			} catch (PointOnEdgeException | NegativeZAxisException
					| PositiveZAxisException e) {
				insideShape += .5;

			}

		}

		double portionOfKnownVolume = insideShape / (double) (coord.size());
		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
	}

	/**
	 * Maintenance method for testing
	 * 'generateUniformDistributionOfSamplePoints' method Method prints a file
	 * containing all the generated sample points
	 */
	public void printUniformDistribution() {
		FileWriter file = null;
		try {
			file = new FileWriter("Distribution");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ArrayList<double[]> coord = generateUniformDistributionOfSamplePoints();

		PrintWriter fileOut = new PrintWriter(file);

		for (int i = 0; i < coord.size(); i++) {
			fileOut.println(String.format("[%.4f, %.4f, %.4f]",
					coord.get(i)[x], coord.get(i)[y], coord.get(i)[z]));
		}

	}

	/**
	 * This method can be used on uniformly distributed sample points. For
	 * random distribution use 'determineIfPointIsInsideShapeRandomDistribution'
	 * Checks if point is between specified bounds for theta. Checks if point is
	 * between specified bounds for phi. Checks if point is within the equation
	 * of an ellipsoid. Throws exceptions for cases where the point is on the
	 * edge of the specified shape.
	 * 
	 * fomulas: theta: thetaStart < thetaSample < thetaEnd phi: phiStart <
	 * phiSample < phiEnd equation of ellipsoid: (x^2/a^2) + (y^2/b^2) +
	 * (z^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	protected boolean determineIfPointIsInsideShapeNonRandomDistribution(
			double[] coord) throws PointOnEdgeException,
			PositiveZAxisException, NegativeZAxisException {

		boolean insideShape = false;

		try {
			if (pointInsideEquationOfEllipsoidRect(coord)
					&& pointBetweenPhiStartAndPhiEndRect(coord)
					&& pointBetweenThetaStartAndThetaEndEvenDist(coord)) {
				insideShape = true;
			}
		} catch (NegativeZAxisException e) {
			if (radianMeasureOffZAxisEnd == Math.PI) {
				throw new NegativeZAxisException(coord);
			}
		} catch (PositiveZAxisException e1) {
			if (radianMeasureOffZAxisStart == 0) {
				throw new PositiveZAxisException(coord);
			}
		}

		return insideShape;
	}

	/**
	 * checks if theta value of sample point is > thetaStart but < thetaEnd. If
	 * theta is == to startTheta or == to endTheta, a PointOnEdgeException is
	 * thrown
	 * 
	 * @param coord
	 * @return
	 * @throws NegativeZAxisException
	 */
	protected boolean pointBetweenThetaStartAndThetaEndEvenDist(double[] coord)
			throws NegativeZAxisException, PositiveZAxisException,
			PointOnEdgeException {
		boolean insideThetaBound = false;

		try {
			if (startRadianTheta <= getThetaValueOfSampleCoordStartTheta(coord)
					&& getThetaValueOfSampleCoordEndTheta(coord) <= endRadianTheta) {

				if (startRadianTheta == getThetaValueOfSampleCoordStartTheta(coord)
						|| getThetaValueOfSampleCoordEndTheta(coord) == endRadianTheta) {
					throw new PointOnEdgeException(coord);
				}

				insideThetaBound = true;
			}
		} catch (PointOnOriginException e) {
			insideThetaBound = true;
			System.out.println(e.getMessage());
		}

		return insideThetaBound;
	}

	/**
	 * Returns an ArrayList<double[]> of 3D sample points of uniform
	 * distribution (equally spaced non random)
	 * 
	 * @return
	 */
	protected ArrayList<double[]> generateUniformDistributionOfSamplePoints() {
		final int samplesOnAxis = 50;

		ArrayList<double[]> coord = new ArrayList<double[]>();
		double xVar = 0;
		for (int xAxis = 0; xAxis < samplesOnAxis; xAxis++) {
			double yVar = 0;
			xVar += a / (double) samplesOnAxis;

			for (int yAxis = 0; yAxis < samplesOnAxis; yAxis++) {
				double zVar = 0;
				yVar += b / (double) samplesOnAxis;

				for (int zAxis = 0; zAxis < samplesOnAxis; zAxis++) {

					zVar += c / (double) samplesOnAxis;

					// first octant
					coord.add(new double[] { xVar, yVar, zVar });
					// second octant
					coord.add(new double[] { xVar * -1, yVar, zVar });
					// third octant
					coord.add(new double[] { xVar * -1, yVar * -1, zVar });
					// fourth octant
					coord.add(new double[] { xVar, yVar * -1, zVar });
					// fifth octant
					coord.add(new double[] { xVar, yVar, zVar * -1 });
					// sixth octant
					coord.add(new double[] { xVar * -1, yVar, zVar * -1 });
					// seventh octant
					coord.add(new double[] { xVar * -1, yVar * -1, zVar * -1 });
					// eighth octant
					coord.add(new double[] { xVar, yVar * -1, zVar * -1 });

				}

			}

		}

		// point on origin
		coord.add(new double[] { 0, 0, 0 });
		// add special conditions
		addZAxisPoints(coord, samplesOnAxis);
		addXYPlanePoints(coord, samplesOnAxis);
		addXAxisPoints(coord, samplesOnAxis);
		addYAxisPoints(coord, samplesOnAxis);
		return coord;

	}

	/**
	 * Helper method for uniform distribution method. Adds points above and
	 * below the origin (z axis) to ArrayList
	 * 
	 * @param coord
	 * @param samplesOnAxis
	 */
	private void addZAxisPoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double zVar = 0;
		zVar += c / (double) samplesOnAxis;

		for (int i = 0; i < samplesOnAxis; i++) {
			zVar += c / (double) samplesOnAxis;

			if (i % 2.0 > 0.0) {
				coord.add(new double[] { 0, 0, zVar });
				coord.add(new double[] { 0, 0, -zVar });
			} else {
				coord.add(new double[] { 0, 0, zVar });
				coord.add(new double[] { 0, 0, -zVar });
			}
		}

	}

	/**
	 * Helper method for uniform distribution method. Adds points contained
	 * within the xy plane to the ArrayList
	 * 
	 * @param coord
	 * @param samplesOnAxis
	 */
	private void addXYPlanePoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double xVar = 0;

		for (int i = 0; i < samplesOnAxis; i++) {

			double yVar = 0;

			for (int index = 0; index < samplesOnAxis; index++) {
				yVar += b / (double) samplesOnAxis;

				// first quadrant
				coord.add(new double[] { xVar, yVar, 0 });
				// second quadrant
				coord.add(new double[] { -xVar, yVar, 0 });
				// third quadrant
				coord.add(new double[] { -xVar, -yVar, 0 });
				// fourth quadrant
				coord.add(new double[] { xVar, -yVar, 0 });
			}

		}
	}

	/**
	 * Helper method for uniform distribution method. Adds points above and
	 * below x axis to ArrayList
	 * 
	 * @param coord
	 * @param samplesOnAxis
	 */
	private void addXAxisPoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double xVar = 0;

		for (int i = 0; i < samplesOnAxis; i++) {
			xVar += a / (double) samplesOnAxis;
			double zVar = 0;

			for (int index = 0; index < samplesOnAxis; index++) {
				zVar += c / (double) samplesOnAxis;

				if (index % 2.0 > 0.0) {
					coord.add(new double[] { xVar, 0, zVar });
					coord.add(new double[] { -xVar, 0, -zVar });
				} else {
					coord.add(new double[] { xVar, 0, zVar });
					coord.add(new double[] { -xVar, 0, -zVar });
				}
			}
		}
	}

	/**
	 * Helper method for uniform distribution method. Adds points above and
	 * below x axis to ArrayList
	 * 
	 * @param coord
	 * @param samplesOnAxis
	 */
	private void addYAxisPoints(ArrayList<double[]> coord, int samplesOnAxis) {
		double yVar = 0;

		for (int i = 0; i < samplesOnAxis; i++) {
			yVar += b / (double) samplesOnAxis;
			double zVar = 0;

			for (int index = 0; index < samplesOnAxis; index++) {
				zVar += c / (double) samplesOnAxis;

				coord.add(new double[] { 0, yVar, zVar });
				coord.add(new double[] { 0, -yVar, -zVar });
			}
		}
	}

}
