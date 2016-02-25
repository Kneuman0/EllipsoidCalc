package biz.personalAcademics.ellipsoidCalc;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleGroup;

public class EllipsoidCalcController {

    @FXML
    private TextField phiStart;
    
    @FXML
    private TextField phiEnd;

    @FXML
    private TextField a;

    @FXML
    private TextField b;

    @FXML
    private TextField c;

    @FXML
    private TextField thetaEnd;

    @FXML
    private Label warningLabel;

    @FXML
    private ToggleGroup angle;

    @FXML
    private RadioButton radiansRadio;

    @FXML
    private Label volumeAnswer;

    @FXML
    private RadioButton degreesRadio;

    @FXML
    private TextField thetaStart;
    
    
    public void initialize(){
    	// Deliberately left blank
    }
    
    public void convertDecimal(ActionEvent e){
    	TextField temp = (TextField) e.getSource();
    	try {
			temp.setText(String.format("%f", EllipsoidCalcUtility.convertToDecimal(temp.getText())));
		} catch (InvalidUserInputException e1) {
			warningLabel.setText(String.format("Non-number %s detected! Fix this to proceed", e1.getMessage()));
		}
    	
    }
    
    public void calculateButton(){
    	String startAngle = thetaStart.getText();
    	String endAngle = thetaEnd.getText();
    	String zAxisAngleEnd = phiEnd.getText();
    	String zAxisAngleStart = phiStart.getText();
    	boolean inDegrees = degreesRadio.isSelected();;
    	double aAxis = 0;
    	double bAxis = 0;
    	double cAxis = 0;
		try {
			aAxis = Double.parseDouble(a.getText());
			bAxis = Double.parseDouble(b.getText());
			cAxis = Double.parseDouble(c.getText());
		} catch (NumberFormatException e1) {
			warningLabel.setText(String.format("Non-number %s detected", e1.getMessage()));
		}
    	try {
			double thetaBegin = EllipsoidCalcUtility.convertThetaToRadians(startAngle, inDegrees);
			double thetaEnd = EllipsoidCalcUtility.convertThetaToRadians(endAngle, inDegrees);
			double phiAngleEnd = EllipsoidCalcUtility.convertThetaToRadians(zAxisAngleEnd, inDegrees);
			double phiAngleStart = EllipsoidCalcUtility.convertThetaToRadians(zAxisAngleStart, inDegrees);
			volumeAnswer.setText(String.format("%f", EllipsoidCalcUtility.getVolume(thetaBegin, thetaEnd, phiAngleEnd, phiAngleStart,
					aAxis, bAxis, cAxis)));
		} catch (InvalidUserInputException e) {
			warningLabel.setText(String.format("Non number %s Detected", e.getMessage()));
		}
    	
    }
    
    
    
}