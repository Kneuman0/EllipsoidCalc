package ellipsoidCalc;

public class TestUtilMethods {

	public static void main(String[] args) throws InvalidUserInputException {
		EllipsoidCalcUtility turd = new EllipsoidCalcUtility();
		
		System.out.println(turd.convertToDecimal("1/2"));
		
		System.out.println(turd.getVolume(Math.PI/4, 3 * Math.PI/4, Math.PI/2, 0, 3, 2, 2));
		
		double three = 3/6; 
		
		System.out.println(three);

	}

}
