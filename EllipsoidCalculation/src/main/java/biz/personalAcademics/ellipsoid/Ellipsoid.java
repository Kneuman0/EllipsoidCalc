package biz.personalAcademics.ellipsoid;

import java.util.Random;

import biz.personalAcademics.ellipsoid.customExceptions.InvalidUserInputException;

public class Ellipsoid extends EllipsoidalShape {

	public static final int MIN_SAMPLE_SIZE = 10_000_000;
	
	// rectangular and cylindrical coordinate indexes
	protected final int x = 0, y = 1, z = 2;
	// spherical and cylindrical coordinate indexes
	protected final int phi = 0;
	protected final int theta = 1;
	protected final int p = 2;
	// cylindrical coordinate indexes
	protected final int r = 0;
	
	protected int sampleSize;
	

	protected Random randomGenerator;

	public Ellipsoid(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisStart, double radianMeasureOffZAxisEnd,
			double a, double b, double c) {

		super(startRadianTheta, endRadianTheta, radianMeasureOffZAxisStart,
				radianMeasureOffZAxisEnd, a, b, c);

		randomGenerator = new Random();
		this.sampleSize = Ellipsoid.MIN_SAMPLE_SIZE;
	}

	/**
	 * Finds the volume of the ellipsoid based off the spherical coordinates the
	 * user passed in. This method is only useful when the values of phi and
	 * theta are increments of pi/2 or is a = b = c. It is entirely inaccurate
	 * otherwise. getEstimatedVolume() should be used in these cases. If volume
	 * is negative (due to incorrectly entered bounds of integration) 0 is
	 * returned.
	 *
	 * @return
	 */
	protected double getExactVolume() {

		double oneThird = 1.0 / 3.0;
		double volume = (-oneThird)
				* (a * b * c * (Math.cos(radianMeasureOffZAxisEnd) - Math
						.cos(radianMeasureOffZAxisStart)))
				* (endRadianTheta - startRadianTheta);

		if (volume < 0) {
			return 0;
		} else {
			return volume;
		}
	}

	/**
	 * Estimates the volume of the ellipsoid using the specified number of
	 * random sample points and the most reliable method developed thus far
	 * 
	 * @param sampleSize
	 * @return
	 */
	public double getEstimatedVolume(int sampleSize)
			throws InvalidUserInputException {
		if (sampleSize < Ellipsoid.MIN_SAMPLE_SIZE) {
			throw new InvalidUserInputException(
					new Integer(sampleSize).toString());
		}

		if (this.executeDefiniteIntegral) {

			return getExactVolume();

		} else {

			if (this.eccentricity >= 1) {
				
				this.sampleSize = sampleSize;
				
				return new EllipsoidSphericalCoords(startRadianTheta,
						endRadianTheta, radianMeasureOffZAxisStart,
						radianMeasureOffZAxisEnd, a, b, c)
						.getEstimatedVolume(sampleSize);

			} else {

				this.sampleSize = sampleSize;
				
				return new EllipsoidRectangularCoords(startRadianTheta,
						endRadianTheta, radianMeasureOffZAxisStart,
						radianMeasureOffZAxisEnd, a, b, c)
						.getEstimatedVolume(sampleSize);

			}

		}
	}

	/**
	 * Estimates the volume of the portion of the ellipsoid using using the most
	 * reliable method developed thus far
	 * 
	 * @return
	 */
	public double getEstimatedVolume() {

		if (this.executeDefiniteIntegral) {

			return getExactVolume();

		} else {

			if (eccentricity >= 1) {

				return new EllipsoidSphericalCoords(startRadianTheta,
						endRadianTheta, radianMeasureOffZAxisStart,
						radianMeasureOffZAxisEnd, a, b, c).getEstimatedVolume();

			} else {

				return new EllipsoidRectangularCoords(startRadianTheta,
						endRadianTheta, radianMeasureOffZAxisStart,
						radianMeasureOffZAxisEnd, a, b, c).getEstimatedVolume();

			}

		}

	}
	
	public double getError(){
		return .0008;
	}

	/**
	 * This method will return the estimated volume of the shape through random
	 * sample points rounded to two decimal places.
	 */
	public String toString() {
		double volume = this.getEstimatedVolume();
		String error;
		double volumeError = getError() * volume;

		if(this.executeDefiniteIntegral){
			error = "0.0  *Definte integral used*";
		}else{			
			error = String.format("%.2f", volumeError);
		}

		return String.format("%.2f +/- %s", volume, error);
	}

	/**
	 * This method will return the estimated volume of the shape through the
	 * specified number of random sample points rounded to two decimal places.
	 */
	public String toString(int sampleSize) {
		double volume = this.getEstimatedVolume(sampleSize);
		String error;
		double volumeError = getError() * volume;

		
		if(this.executeDefiniteIntegral){
			error = "0.0  *Definte integral used*";
		}else{			
			error = String.format("%.2f", volumeError);
		}

		return String.format("%.2f +/- %s", volume, error);
	}

	/**
	 * Returns a positive 1 or -1 randomly
	 * 
	 * @return
	 */
	protected int getRandomNegation() {
		boolean negativeSign = randomGenerator.nextBoolean();

		if (negativeSign) {
			return 1;
		} else {
			return -1;
		}
	}
	
	protected int getSampleSize(){
		return sampleSize;
	}
}
