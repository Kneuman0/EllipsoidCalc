package biz.personalAcademics.ellipsoid;

import biz.personalAcademics.ellipsoid.customExceptions.InvalidUserInputException;
import static java.lang.Math.PI;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static java.lang.Math.sin;
import static java.lang.Math.cos;


public class EllipsoidCylindricalCoords extends Ellipsoid{

		
	public EllipsoidCylindricalCoords(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisStart, double radianMeasureOffZAxisEnd,
			double a, double b, double c) {
		
		super(startRadianTheta, endRadianTheta, radianMeasureOffZAxisStart, radianMeasureOffZAxisEnd,
			a, b, c);
	}
	
	/**
	 * Estimates the volume of the ellipsoid using 5,000,000 random sample
	 * points. The algorithm uses cylindrical coordinates to generate random
	 * sample points. The sample points will be contained inside the known
	 * volume of a cylinder defined by:
	 * 
	 * (pi)(r^2)(h) = (pi) * (a^2) * (2 * c).
	 * 
	 * a MUST be the major axis in the xy plane and is associated with the
	 * positive x axis.
	 * 
	 * c MUST be associated with the z axis.
	 * 
	 * @param sampleSize
	 * @return
	 */
	public double getEstimatedVolume() {
		int insideShape = 0;
		final int sampleSize = 5_000_000;

		for (int i = 0; i < sampleSize; i++) {
			if (pointInsideShapeCyl(generate3DSamplePointCyl())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system V = (pi)(h)(r^2)
		return portionOfKnownVolume * PI * pow(this.getA(), 2) * 2
				* this.getC();
	}

	/**
	 * Estimates the volume of the ellipsoid using the specified number of
	 * random sample points. The algorithm uses cylindrical coordinates to
	 * generate random sample points. The sample points will be contained inside
	 * the known volume of a cylinder defined by:
	 * 
	 * (pi)(r^2)(h) = (pi) * (a^2) * (2 * c).
	 * 
	 * a MUST be the major axis in the xy plane and is associated with the
	 * positive x axis.
	 * 
	 * c MUST be associated with the z axis.
	 * 
	 * @param sampleSize
	 * @return
	 */
	public double getEstimatedVolume(int sampleSize) throws InvalidUserInputException {
		if (sampleSize < Ellipsoid.MIN_SAMPLE_SIZE){
			throw new InvalidUserInputException(new Integer(sampleSize).toString());
		}
		
		int insideShape = 0;

		for (int i = 0; i < sampleSize; i++) {
			if (pointInsideShapeCyl(generate3DSamplePointCyl())) {
				insideShape++;
			}

		}

		double portionOfKnownVolume = insideShape / (double) sampleSize;

		// portion of volume of 3D system V = (pi)(h)(r^2)
		return portionOfKnownVolume * PI * pow(this.getA(), 2) * 2
				* this.getC();
	}
	
	/**
	 * Checks if point is inside equation of ellipsoid, between the values of
	 * phi and inside the equation of the ellipse in the xy plane
	 * 
	 * @param coord
	 * @return
	 */
	private boolean pointInsideShapeCyl(double[] coord) {

		return checkIfThetaIsBetweenStartAndEndCyl(coord)
				&& pointBetweenPhiValuesCyl(coord)
				&& pointInsideEquationOfEllipsoidCyl(coord)
				&& checkLengthOfR(coord);
	}

	/**
	 * calculates the longest possible value of r with the given values of
	 * theta. If the random value of r is longer than the calculated value, the
	 * point is outside the shape.
	 * 
	 * @param coord
	 * @return
	 */
	private boolean checkLengthOfR(double[] coord) {
		double numerator = this.getA() * this.getB();

		// ( b^2 cos^2(theta) + a^2 sin^2(theta) )^(1/2)
		double denominator = sqrt(pow(
				this.getB() * cos(coord[theta]), 2)
				+ pow(a * sin(coord[theta]), 2));
		double rValue = numerator / denominator;

		return coord[r] <= rValue;

	}

	/**
	 * calculates: r^2*cos^2(theta)/a^2 + r^2*sin^2(theta)/b^2 + z^2/c^2
	 * 
	 * then checks if that value is <= to 1. If not, the point is outside of the
	 * shape
	 * 
	 * @param coord
	 * @return
	 */
	private boolean pointInsideEquationOfEllipsoidCyl(double[] coord) {
		double xValue = pow((coord[r] * cos(coord[theta])) / a,
				2);
		double yValue = pow((coord[r] * sin(coord[theta])) / b,
				2);
		double zValue = pow(coord[z] / c, 2);

		double pointValue = xValue + yValue + zValue;

		return pointValue <= 1;
	}

	/**
	 * Checks if the the random z value is below the value of phiStart and above
	 * the values of phiEnd. If not, the point is outside the shape.
	 * 
	 * @param coord
	 * @return
	 */
	private boolean pointBetweenPhiValuesCyl(double[] coord) {

		return checkThatZMoreThanPhiLineEnd(coord)
				&& checkThatZLessThanPhiLineStart(coord);
	}

	private double cot(double theta) {
		return (cos(theta) / sin(theta));
	}

	/**
	 * Calculates: z = cot(phiEnd) * r
	 * 
	 * Then checks if the random value of z is greater than the calculated
	 * value. If not, the point is outside of the shape.
	 * 
	 * This method assumes phiEnd is never 0
	 * 
	 * @param coord
	 * @return
	 */
	private boolean checkThatZMoreThanPhiLineEnd(double[] coord) {
		if (radianMeasureOffZAxisEnd == PI) {

			return true;
		} else {

			double phiValue = cot(radianMeasureOffZAxisEnd) * coord[r];

			if (phiValue <= coord[z]) {

				return true;

			} else {

				return false;
			}
		}
	}

	/**
	 * Calculates: z = cot(phiEnd) * r
	 * 
	 * Then checks if the random value of z is less than the calculated value.
	 * If not, the point is outside of the shape.
	 * 
	 * This method assumes phiStart is never pi
	 * 
	 * @param coord
	 * @return
	 */
	private boolean checkThatZLessThanPhiLineStart(double[] coord) {
		if (radianMeasureOffZAxisStart == 0) {

			return true;
		} else {

			double phiValue = (cos(radianMeasureOffZAxisStart) / 
					sin(radianMeasureOffZAxisStart)) * coord[r];

			if (coord[z] <= phiValue) {

				return true;

			} else {

				return false;
			}
		}
	}

	/**
	 * Checks thetaStart <= theta <= thetaEnd
	 * 
	 * If not, the point is outside the shape.
	 * 
	 * @param coord
	 * @return
	 */
	private boolean checkIfThetaIsBetweenStartAndEndCyl(double[] coord) {
		boolean insideThetaBound = false;

		if (startRadianTheta <= coord[theta] && coord[theta] <= endRadianTheta) {

			insideThetaBound = true;
		}

		return insideThetaBound;
	}

	/**
	 * Generates a random point inside of a cylinder with uniform distribution
	 * using cylindrical coordinates
	 * 
	 * index 0 = r, index 1 = theta, index 2 = z
	 * 
	 * @return
	 */
	private double[] generate3DSamplePointCyl() {
		double[] samplePoint = new double[3];

		// r will range from 0 to the longest axis
		samplePoint[r] = this.a * sqrt(randomGenerator.nextDouble());

		// theta will range 0 to 2pi
		samplePoint[theta] = 2 * PI * randomGenerator.nextDouble();

		// z will range from 0 to the c axis both above and below the xy plane
		samplePoint[z] = this.c * randomGenerator.nextDouble()
				* getRandomNegation();

		return samplePoint;
	}

}
