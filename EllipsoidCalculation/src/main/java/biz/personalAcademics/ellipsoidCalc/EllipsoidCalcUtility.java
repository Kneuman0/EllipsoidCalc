package biz.personalAcademics.ellipsoidCalc;

public class EllipsoidCalcUtility {
	
	/**
	 * Finds the volume of the ellipsoid based off the spherical coordinates the user passed in
	 * @param startRadianTheta
	 * @param endRadianTheta
	 * @param radianMeasureOffZAxisEnd
	 * @param radianMeasureOffZAxisStart
	 * @param a
	 * @param b
	 * @param c
	 * @return
	 */
	public static double getVolume(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisEnd, double radianMeasureOffZAxisStart,
			double a, double b, double c) {
		
		double oneThird = 1 / (double)3;
		double volume = (-oneThird)*(a*b*c*(Math.cos(radianMeasureOffZAxisEnd) - 
				Math.cos(radianMeasureOffZAxisStart)))*(endRadianTheta - startRadianTheta);
		return volume;
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
	 * If user inputs angle in degrees, it is converted and returned in radians.
	 * If radians were chosen, pi is multiplied in 
	 * @param degrees
	 * @param degreeIsSelected
	 * @return
	 * @throws InvalidUserInputException
	 */
	public static double convertThetaToRadians(String degrees, boolean degreeIsSelected)
			throws InvalidUserInputException {
		double testDegree;

		try {
			testDegree = Double.parseDouble(degrees);
		} catch (NumberFormatException e1) {
			throw new InvalidUserInputException(degrees);
		}

		if (degreeIsSelected) {
			return (testDegree / 180) * (Math.PI);
		} else {
			return testDegree * Math.PI;
		}
	}

}
