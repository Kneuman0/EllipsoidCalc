package ellipsoidCalc;

public class InvalidUserInputException extends RuntimeException {
     public InvalidUserInputException(String invalidInput){
    	 super(invalidInput);
     }
}
