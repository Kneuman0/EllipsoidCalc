package biz.personalAcademics.ellipsoidCalc;

@SuppressWarnings("serial")
public class InvalidUserInputException extends RuntimeException {
     public InvalidUserInputException(String invalidInput){
    	 super(invalidInput);
     }
}
