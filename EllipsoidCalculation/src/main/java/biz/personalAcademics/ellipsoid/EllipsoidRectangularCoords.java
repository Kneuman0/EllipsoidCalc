package biz.personalAcademics.ellipsoid;

import biz.personalAcademics.ellipsoid.customExceptions.*;

public class EllipsoidRectangularCoords extends Ellipsoid{
	
	public EllipsoidRectangularCoords(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisStart, double radianMeasureOffZAxisEnd,
			double a, double b, double c) {
		
		super(startRadianTheta, endRadianTheta, radianMeasureOffZAxisStart, radianMeasureOffZAxisEnd,
			a, b, c);
	}
	
	/**
	 * creates 10 million sample points inside a rectangular prism using a
	 * rectangular coordinate system the rectangular prism's volume is defined
	 * as 2a * 2b * 2c where a b and c represent the axes of the ellipsoid. This
	 * method inscribes an ellipsoid inside the rectangular prism to compute the
	 * estimated volume. The error has empirically been found to be as much as
	 * +/- .1
	 * 
	 * @return
	 */
	public double getEstimatedVolume() {
		int insideShape = 0;
		final int sampleSize = 10_000_000;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShapeRandomDistributionRect(generateRandom3DSamplePointRect())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
	}

	/**
	 * creates a specified number of sample points inside a rectangular prism
	 * using a rectangular coordinate system the rectangular prism's volume is
	 * defined as 2a * 2b * 2c where a b and c represent the axes of the
	 * ellipsoid. This method inscribes an ellipsoid inside the rectangular
	 * prism to compute the estimated volume. The error has empirically been
	 * found to be as much as +/- .1
	 * 
	 * @return
	 */
	public double getEstimatedVolume(int sampleSize) throws InvalidUserInputException{
		if (sampleSize < Ellipsoid.MIN_SAMPLE_SIZE){
			throw new InvalidUserInputException(new Integer(sampleSize).toString());
		}
		
		int insideShape = 0;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShapeRandomDistributionRect(generateRandom3DSamplePointRect())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system
		return portionOfKnownVolume * (2 * a) * (2 * b) * (2 * c);
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
	protected boolean determineIfPointIsInsideShapeRandomDistributionRect(
			double[] coord) {

		return pointInsideEquationOfEllipsoidRect(coord)
				&& pointBetweenPhiStartAndPhiEndRect(coord)
				&& pointBetweenThetaStartAndThetaEndRandomDistRect(coord);
	}

	/**
	 * checks if sample point is inside the equation of an ellipsoid using the
	 * boolean expression: (x^2/a^2) + (y^2/b^2) + (z^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	boolean pointInsideEquationOfEllipsoidRect(double[] coord) {
		// point is inside shape if x^2/a^2 + y^2/b^2 + z^2/c^2 <= 1 then it is
		// inside the ellipsoid
		double valueOfPointInEllipsoidEquation = (coord[x] * coord[x])
				/ (a * a) + (coord[y] * coord[y]) / (b * b)
				+ (coord[z] * coord[z]) / (c * c);

		return valueOfPointInEllipsoidEquation <= 1;
	}

	/**
	 * Checks if the phi value of sample point is >= phiStart but <= phiEnd
	 * 
	 * @param coord
	 * @return
	 */
	protected boolean pointBetweenPhiStartAndPhiEndRect(double[] coord)
			throws PointOnEdgeException {
		boolean insidePhiBound = false;

		try {
			if (radianMeasureOffZAxisStart <= getPhiValueOfSampleCoord(coord)
					&& getPhiValueOfSampleCoord(coord) <= radianMeasureOffZAxisEnd) {

				if (radianMeasureOffZAxisStart == getPhiValueOfSampleCoord(coord)
						|| getPhiValueOfSampleCoord(coord) == radianMeasureOffZAxisEnd) {

					throw new PointOnEdgeException(coord);
				}

				insidePhiBound = true;

			}
		} catch (PointOnOriginException e) {
			insidePhiBound = true;
			System.out.println(e.getMessage());
		}
		return insidePhiBound;
	}

	/**
	 * checks if theta value of sample point is > thetaStart but < thetaEnd.
	 * Randomly decides whether or not to include special cases.
	 * 
	 * @param coord
	 * @return
	 * @throws NegativeZAxisException
	 */
	protected boolean pointBetweenThetaStartAndThetaEndRandomDistRect(
			double[] coord) {
		boolean insideThetaBound = false;

		try {
			if (startRadianTheta <= getThetaValueOfSampleCoordStartTheta(coord)
					&& getThetaValueOfSampleCoordEndTheta(coord) <= endRadianTheta) {

				insideThetaBound = true;
			}
		} catch (NegativeArraySizeException | PointOnOriginException
				| PositiveZAxisException e) {
			System.out.println(e.getMessage());
			// Randomly determines whether or not to include point in shape
			insideThetaBound = randomGenerator.nextBoolean();
		}

		return insideThetaBound;
	}

	/**
	 * Calculates the value of phi using the formula: arcCos(z/sqrt(x^2 + y^2 +
	 * z^2)). If z is negative, pi/2 is added to the angle. If z = 0 but x != 0
	 * or y != 0 then phi is pi/2. If x == 0 and y == 0 and z != 0, then phi is
	 * either 0 or pi. If x == 0, y == 0, z == 0, then PointOnOrigin exception
	 * is thrown.
	 * 
	 * @param coord
	 * @return
	 * @throws PointOnOriginException
	 */
	protected double getPhiValueOfSampleCoord(double[] coord)
			throws PointOnOriginException {
		// above xy plane (positive z value)
		if (coord[z] > 0 && coord[x] != 0 && coord[y] != 0) {

			// arcCos(z/sqrt(x^2 + y^2 + z^2))
			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// below xy plane (negative z value)
		} else if (coord[z] < 0 && coord[x] != 0 && coord[y] != 0) {

			// arcCos(z/sqrt(x^2 + y^2 + z^2)) + pi
			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])))
			// + Math.PI
			;

			// on positive z axis = 0
		} else if (coord[z] > 0 && coord[x] == 0 && coord[y] == 0) {
			return 0;

			// on negative z axis = pi
		} else if (coord[z] < 0 && coord[x] == 0 && coord[y] == 0) {
			return Math.PI;

			// in zy plane above xy plane
		} else if (coord[z] > 0 && coord[x] == 0 && coord[y] != 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// in zy plane below xy plane
		} else if (coord[z] < 0 && coord[x] == 0 && coord[y] != 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// in zx plane above xy plane
		} else if (coord[z] > 0 && coord[x] != 0 && coord[y] == 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// in zx plane below xy plane
		} else if (coord[z] < 0 && coord[x] != 0 && coord[y] == 0) {

			return Math.acos(coord[z]
					/ (Math.sqrt(coord[x] * coord[x] + coord[y] * coord[y]
							+ coord[z] * coord[z])));

			// in xy plane
		} else if (coord[z] == 0) {

			if (coord[x] != 0 || coord[y] != 0) {
				return Math.PI / 2;

			} else {

				throw new PointOnOriginException(coord);

			}

			// on origin
		} else {

			throw new PointOnOriginException(coord);

		}

		// do case where phi is below x y plane

	}

	/**
	 * Returns the values of theta in the x y plane using the formula
	 * arcTan(y/x). ***If point is on positive x axis then 0 is returned.*** If
	 * point is in first quadrant then the angle is returned. If point is in
	 * second quadrant, then pi/2 is added to angle. If point is in third
	 * quadrant, then pi/2 is added to angle. If point is in fourth quadrant
	 * then 2pi is added to the angle. If x = 0, y = 0 and z > 0, then
	 * PositiveZAxisException is thrown. If x = 0, y = 0 and z < 0 then
	 * NegativeArraySizeException is thrown. If x = 0, y = 0 and z = 0, then
	 * PointOnOriginException is thrown.
	 * 
	 * @param coord
	 * @return
	 * @throws PositiveZAxisException
	 * @throws NegativeArraySizeException
	 * @throws PointOnOriginException
	 */
	protected double getThetaValueOfSampleCoordStartTheta(double[] coord)
			throws PositiveZAxisException, NegativeArraySizeException,
			PointOnOriginException {
		// first quadrant
		if (coord[x] > 0 && coord[y] > 0) {
			return Math.atan(coord[y] / coord[x]);

			// Second quadrant
		} else if (coord[x] < 0 && coord[y] > 0) {
			return Math.atan(coord[y] / coord[x]) + Math.PI;

			// Third Quadrant
		} else if (coord[x] < 0 && coord[y] < 0) {
			return Math.atan(coord[y] / coord[x]) + Math.PI;

			// Fourth Quadrant
		} else if (coord[x] > 0 && coord[y] < 0) {
			return Math.atan(coord[y] / coord[x]) + 2 * Math.PI;

			// on positive y axis
		} else if (coord[x] == 0 && coord[y] > 0) {
			return Math.PI / 2;

			// on negative y axis
		} else if (coord[x] == 0 && coord[y] < 0) {
			return (3 / 2) * Math.PI;

			// on positive x axis (choosing 0 for radian measure because this is
			// the start value)
		} else if (coord[x] > 0 && coord[y] == 0) {
			return 0;

			// on negative x axis
		} else if (coord[x] < 0 && coord[y] == 0) {
			return Math.PI;

			// on negative z axis
		} else if (coord[x] == 0 && coord[y] == 0 && coord[z] < 0) {

			throw new NegativeZAxisException(coord);

		} else if (coord[x] == 0 && coord[y] == 0 && coord[z] > 0) {
			throw new PositiveZAxisException(coord);

			// on positive z axis
		} else {
			throw new PointOnOriginException(coord);
		}
	}

	/**
	 * Returns the values of theta in the x y plane using the formula
	 * arcTan(y/x). ***If point is on positive x axis then 2pi is returned.***
	 * If point is in first quadrant then the angle is returned. If point is in
	 * second quadrant, then pi/2 is added to angle. If point is in third
	 * quadrant, then pi/2 is added to angle. If point is in fourth quadrant
	 * then 2pi is added to the angle. If x = 0, y = 0 and z > 0, then
	 * PositiveZAxisException is thrown. If x = 0, y = 0 and z < 0 then
	 * NegativeArraySizeException is thrown. If x = 0, y = 0 and z = 0, then
	 * PointOnOriginException is thrown.
	 * 
	 * @param coord
	 * @return
	 * @throws PositiveZAxisException
	 * @throws NegativeArraySizeException
	 * @throws PointOnOriginException
	 */
	protected double getThetaValueOfSampleCoordEndTheta(double[] coord)
			throws PositiveZAxisException, NegativeArraySizeException,
			PointOnOriginException {
		// first quadrant
		if (coord[x] > 0 && coord[y] > 0) {
			return Math.atan(coord[y] / coord[x]);

			// Second quadrant
		} else if (coord[x] < 0 && coord[y] > 0) {
			return Math.atan(coord[y] / coord[x]) + Math.PI;

			// Third Quadrant
		} else if (coord[x] < 0 && coord[y] < 0) {
			return Math.atan(coord[y] / coord[x]) + Math.PI;

			// Fourth Quadrant
		} else if (coord[x] > 0 && coord[y] < 0) {
			return Math.atan(coord[y] / coord[x]) + 2 * Math.PI;

			// on positive y axis
		} else if (coord[x] == 0 && coord[y] > 0) {
			return Math.PI / 2;

			// on negative y axis
		} else if (coord[x] == 0 && coord[y] < 0) {
			return (3 / 2) * Math.PI;

			// on positive x axis (choosing 0 for radian measure because this is
			// the start value)
		} else if (coord[x] > 0 && coord[y] == 0) {
			return 2 * Math.PI;

			// on negative x axis
		} else if (coord[x] < 0 && coord[y] == 0) {
			return Math.PI;

			// on negative z axis
		} else if (coord[x] < 0 && coord[y] == 0 && coord[z] < 0) {

			throw new NegativeZAxisException(coord);

			// on positive z axis
		} else if (coord[x] < 0 && coord[y] == 0 && coord[z] > 0) {

			throw new PositiveZAxisException(coord);

			// on origin
		} else {
			throw new PointOnOriginException(coord);
		}
	}



	/**
	 * generates a 3 dimensional sample point with uniform distribution that
	 * lays somewhere within the three dimensional coordinate system where each
	 * axis ranges from x = +/- a, y = +/- b and z = +/- c
	 * 
	 * index 0 = x, index 1 = y, index 2 = z
	 * 
	 * @return
	 */
	protected double[] generateRandom3DSamplePointRect() {
		double[] samplePoint = new double[3];

		samplePoint[x] = a * randomGenerator.nextDouble() * getRandomNegation();
		samplePoint[y] = b * randomGenerator.nextDouble() * getRandomNegation();
		samplePoint[z] = c * randomGenerator.nextDouble() * getRandomNegation();

		return samplePoint;
	}

}
