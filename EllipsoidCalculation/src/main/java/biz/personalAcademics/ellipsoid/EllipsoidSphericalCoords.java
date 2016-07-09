package biz.personalAcademics.ellipsoid;

import biz.personalAcademics.ellipsoid.customExceptions.InvalidUserInputException;
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.cos;
import static java.lang.Math.acos;
import static java.lang.Math.cbrt;

public class EllipsoidSphericalCoords extends Ellipsoid{
	
	public EllipsoidSphericalCoords(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisStart, double radianMeasureOffZAxisEnd,
			double a, double b, double c) {
		
		super(startRadianTheta, endRadianTheta, radianMeasureOffZAxisStart, radianMeasureOffZAxisEnd,
			a, b, c);
	}
	
	/**
	 * creates a specified number of sample points inside a sphere using a
	 * spherical coordinate system the sphere's volume is defined as
	 * (4/3)pi(largestAxis)^3. This method inscribes an ellipsoid inside the
	 * sphere to compute the estimated volume. The error has empirically been
	 * found to be as much as +/- .03
	 * 
	 * @return
	 */
	public double getEstimatedVolume() {
		int insideShape = 0;
		final int sampleSize = 10_000_000;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShapeRandomDistributionSphere(generateRandom3DSamplePointSphere())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system
		return portionOfKnownVolume * (4.0 / 3.0) * PI
				* pow(getLongestAxis(), 3);
	}

	/**
	 * creates 30 million sample points inside a sphere using a spherical
	 * coordinate system the sphere's volume is defined as
	 * (4/3)pi(largestAxis)^3. This method inscribes an ellipsoid inside the
	 * sphere to compute the estimated volume. The error has empirically been
	 * found to be as much as +/- .03
	 * 
	 * @return
	 */
	public double getEstimatedVolume(int sampleSize) throws InvalidUserInputException{
		if (sampleSize < Ellipsoid.MIN_SAMPLE_SIZE){
			throw new InvalidUserInputException(new Integer(sampleSize).toString());
		}
		
		int insideShape = 0;

		for (int i = 0; i < sampleSize; i++) {
			if (determineIfPointIsInsideShapeRandomDistributionSphere(generateRandom3DSamplePointSphere())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system
		return portionOfKnownVolume * (4.0 / 3.0) * PI
				* pow(getLongestAxis(), 3);
	}

	
	/**
	 * This method can only be used for Spherical coordinates. Checks if point
	 * is between specified bounds for phi. Checks if point is within the
	 * equation of an ellipsoid.
	 * 
	 * fomulas: theta: thetaStart < thetaSample < thetaEnd phi: phiStart <
	 * phiSample < phiEnd equation of ellipsoid: (x^2/a^2) + (y^2/b^2) +
	 * (z^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	private boolean determineIfPointIsInsideShapeRandomDistributionSphere(
			double[] coord) {

		return pointInsideEquationOfEllipsoidSphere(coord)
				&& pointBetweenPhiStartAndPhiEndSphere(coord)
				&& pointBetweenThetaStartAndThetaEndRandomDistSphere(coord);
	}

	/**
	 * Checks if the phi value of sample point is >= phiStart but <= phiEnd
	 * 
	 * @param coord
	 * @return
	 */
	private boolean pointBetweenPhiStartAndPhiEndSphere(double[] coord) {

		return radianMeasureOffZAxisStart <= coord[super.phi]
				&& coord[super.phi] <= radianMeasureOffZAxisEnd;
	}

	/**
	 * checks if theta value of sample point is > thetaStart but < thetaEnd.
	 * Randomly decides whether or not to include special cases.
	 * 
	 * @param coord
	 * @return
	 * @throws NegativeZAxisException
	 */
	protected boolean pointBetweenThetaStartAndThetaEndRandomDistSphere(
			double[] coord) {

		return startRadianTheta <= coord[super.theta]
				&& coord[super.theta] <= endRadianTheta;
	}

	/**
	 * checks if sample point is inside the equation of an ellipsoid using the
	 * boolean expression: ((psin(phi)cos(theta))^2/a^2) +
	 * ((psin(phi)cos(theta))^2/b^2) + ((pcos(phi))^2/c^2) <= 1
	 * 
	 * @param coord
	 * @return
	 */
	protected boolean pointInsideEquationOfEllipsoidSphere(double[] coord) {
		double xSquared = pow(
				(coord[super.p] * sin(coord[super.phi]) * cos(coord[super.theta])), 2);
		double ySquared = pow(
				(coord[super.p] * sin(coord[super.phi]) * sin(coord[super.theta])), 2);
		double zSquared = pow((coord[super.p] * cos(coord[super.phi])), 2);

		/*
		 * 
		 * point is inside shape if ((p * sin(phi)cos(theta))^2/a^2) + ((p *
		 * sin(phi)cos(theta))^2/b^2) + ((p * cos(phi))^2/c^2) <= 1 then it is
		 * inside the ellipsoid
		 */
		double valueOfPointInEllipsoidEquation = (xSquared) / (a * a)
				+ (ySquared) / (b * b) + (zSquared) / (c * c);

		return valueOfPointInEllipsoidEquation <= 1;
	}

	/**
	 * generates a random point in 3 dimensions with uniform distribution using
	 * spherical coordinates. the point lies somewhere inside the sphere
	 * (4/3)pi(r^2) where r is the longest axis.
	 * 
	 * index 0 = phi, index 1 = theta, index 2 = p
	 * 
	 * @return
	 */
	protected double[] generateRandom3DSamplePointSphere() {
		double[] samplePoint = new double[3];

		// phi ranges from 0 to pi
		samplePoint[super.phi] = acos(2 * randomGenerator.nextDouble() - 1);

		// theta ranges from 0 to 2pi
		samplePoint[super.theta] = 2 * PI * randomGenerator.nextDouble();

		// p ranges from 0 to the longest axis
		samplePoint[super.p] = getLongestAxis()
				* cbrt(randomGenerator.nextDouble());

		return samplePoint;
	}
}
