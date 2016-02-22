package ellipsoidCalc;

public class InvalidUserInputException extends Exception {
     public InvalidUserInputException(String invalidInput){
    	 super(invalidInput);
     }
}
