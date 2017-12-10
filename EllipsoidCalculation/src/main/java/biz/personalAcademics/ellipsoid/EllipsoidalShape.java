package biz.personalAcademics.ellipsoid;

import java.util.Arrays;
import static java.lang.Math.PI;

import biz.personalAcademics.ellipsoid.customExceptions.InvalidUserInputException;

public abstract class EllipsoidalShape {
	protected double startRadianTheta, endRadianTheta, radianMeasureOffZAxisEnd,
			radianMeasureOffZAxisStart, a, b, c;
	
	protected boolean executeDefiniteIntegral;

	protected double[] sortedAxes;
	
	protected double radianSum, eccentricity;

	public EllipsoidalShape(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisStart, double radianMeasureOffZAxisEnd,
			double a, double b, double c) {

		this.startRadianTheta = startRadianTheta;
		this.endRadianTheta = endRadianTheta;
		this.radianMeasureOffZAxisEnd = radianMeasureOffZAxisEnd;
		this.radianMeasureOffZAxisStart = radianMeasureOffZAxisStart;
		this.a = a;
		this.b = b;
		this.c = c;
		
		updateFields();
	}
	
	/**
	 * If a setter changes a field, this method will execute to maintain proper
	 * execution pattern for volume calculation
	 */
	private void updateFields() {
		sortedAxes = new double[] { this.a, this.b, this.c };
		Arrays.sort(sortedAxes);
		
		this.eccentricity = (this.c / this.a) + (this.b / this.a);
		
		this.radianSum = this.radianMeasureOffZAxisEnd
				+ this.radianMeasureOffZAxisStart + this.startRadianTheta
				+ this.endRadianTheta;

		this.executeDefiniteIntegral = this.radianSum % (PI / 2) == 0
				|| this.eccentricity == 2;
	}
	

	/**
	 * @return the startRadianTheta
	 */
	public double getStartRadianTheta() {
		return startRadianTheta;
	}

	/**
	 * @param startRadianTheta
	 *            the startRadianTheta to set
	 */
	public void setStartRadianTheta(double startRadianTheta) {
		this.startRadianTheta = startRadianTheta;
		updateFields();
	}

	/**
	 * @return the endRadianTheta
	 */
	public double getEndRadianTheta() {
		return endRadianTheta;
	}

	/**
	 * @param endRadianTheta
	 *            the endRadianTheta to set
	 */
	public void setEndRadianTheta(double endRadianTheta) {
		this.endRadianTheta = endRadianTheta;
		updateFields();
	}

	/**
	 * @return the radianMeasureOffZAxisEnd
	 */
	public double getRadianMeasureOffZAxisEnd() {
		return radianMeasureOffZAxisEnd;
	}

	/**
	 * @param radianMeasureOffZAxisEnd
	 *            the radianMeasureOffZAxisEnd to set
	 */
	public void setRadianMeasureOffZAxisEnd(double radianMeasureOffZAxisEnd) {
		this.radianMeasureOffZAxisEnd = radianMeasureOffZAxisEnd;
		updateFields();
	}

	/**
	 * @return the radianMeasureOffZAxisStart
	 */
	public double getRadianMeasureOffZAxisStart() {
		return radianMeasureOffZAxisStart;
	}

	/**
	 * @param radianMeasureOffZAxisStart
	 *            the radianMeasureOffZAxisStart to set
	 */
	public void setRadianMeasureOffZAxisStart(double radianMeasureOffZAxisStart) {
		this.radianMeasureOffZAxisStart = radianMeasureOffZAxisStart;
		updateFields();
	}
	

	/**
	 * returns the eccentricity of the whole ellipsoid. Formula for
	 * eccentricity: c/a + b/a
	 * 
	 * IF eccentricity = 2, the ellipsoid is a sphere
	 * 
	 * @return
	 */
	public double getTotalEccentricity() {
		return this.eccentricity;
	}

	/**
	 * If definite integral will execute, method will return true otherwise,
	 * method will return false indicating the Monte Carlo integration will be
	 * used.
	 * 
	 * @return
	 */
	public boolean determineIfDefiniteIntegralWillExecute() {
		return executeDefiniteIntegral;
	}

	/**
	 * @return the a
	 */
	public double getA() {
		return a;
	}

	/**
	 * @param a
	 *            the a to set
	 */
	public void setA(double a) {
		this.a = a;
		updateFields();
	}

	/**
	 * @return the b
	 */
	public double getB() {
		return b;
	}

	/**
	 * @param b
	 *            the b to set
	 */
	public void setB(double b) {
		this.b = b;
		updateFields();
	}

	/**
	 * @return the c
	 */
	public double getC() {
		return c;
	}

	/**
	 * @param c
	 *            the c to set
	 */
	public void setC(double c) {
		this.c = c;
		updateFields();
	}

	/**
	 * Converts value in text to a decimal
	 * 
	 * @param input
	 * @return
	 * @throws InvalidUserInputException
	 */
	public static double convertToDecimal(String input)
			throws InvalidUserInputException {
		String[] fraction = input.split("/");
		// Returns user input if array is not a fraction
		if (fraction.length == 1) {
			try {
				return Double.parseDouble(input);
			} catch (NumberFormatException e) {
				throw new InvalidUserInputException(input);
			}
			// converts user input into a decimal if a fraction or throws
			// invalid input exception
		} else {
			int[] numbers = new int[2];
			for (int i = 0; i < fraction.length; i++) {
				try {
					numbers[i] = Integer.parseInt(fraction[i]);
				} catch (NumberFormatException | ArrayIndexOutOfBoundsException e) {
					throw new InvalidUserInputException(input);
				}
			}
			return numbers[0] / (double) numbers[1];
		}

	}

	/**
	 * Returns the longest axis
	 * @return
	 */
	protected double getLongestAxis() {
		return sortedAxes[2];
	}
	
	/**
	 * Returns an array of the given sorted axes smallest to largest
	 * @return
	 */
	public double[] getSortedAxes(){
		return sortedAxes;
	}

	/**
	 * If user inputs angle in degrees, it is converted and returned in radians.
	 * If radians were chosen, pi is multiplied in
	 * 
	 * @param degrees
	 * @param degreeIsSelected
	 * @return
	 * @throws InvalidUserInputException
	 */
	public static double convertThetaToRadians(String degrees,
			boolean degreeIsSelected) throws InvalidUserInputException {
		double testDegree;

		try {
			testDegree = Double.parseDouble(degrees);
		} catch (NumberFormatException e1) {

			// if in fraction form, it is converted to decimal, otherwise
			// invalidUserInputException is thrown
			testDegree = Ellipsoid
					.convertToDecimal(degrees.replaceAll(" ", ""));
		}

		if (degreeIsSelected) {
			return (testDegree / 180) * (Math.PI);
		} else {
			return testDegree * Math.PI;
		}
	}
}
