package ellipsoidCalc;

public class EllipsoidCalcUtility {
	public double getVolume(double startRadianTheta, double endRadianTheta,
			double radianMeasureOffZAxisEnd, double radianMeasureOffZAxisStart,
			double a, double b, double c) {
		
		double oneThird = 1 / (double)3;
		double volume = (-oneThird * (Math.cos(radianMeasureOffZAxisEnd)
				* endRadianTheta * a * b * c)) + (oneThird * (endRadianTheta * a * b * c) * Math.cos(radianMeasureOffZAxisStart))
						+ (oneThird * (Math.cos(radianMeasureOffZAxisEnd) * startRadianTheta * a * b * c))
						- (oneThird * startRadianTheta * a * b * c * Math.cos(radianMeasureOffZAxisStart));
		return volume * -1;
	}

	/**
	 * Converts value in text to a decimal
	 * 
	 * @param input
	 * @return
	 * @throws InvalidUserInputException
	 */
	public double convertToDecimal(String input)
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

	public double convertThetaToRadians(String degrees, boolean degreeIsSelected)
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
