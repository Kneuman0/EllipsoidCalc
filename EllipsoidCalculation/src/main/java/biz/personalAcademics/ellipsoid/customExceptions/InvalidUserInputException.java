package biz.personalAcademics.ellipsoid.customExceptions;

@SuppressWarnings("serial")
public class InvalidUserInputException extends RuntimeException {
     public InvalidUserInputException(String invalidInput){
    	 super(invalidInput);
     }
}
