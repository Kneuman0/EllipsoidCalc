package biz.personalAcademics.ellipsoidCalc;

import org.omg.PortableInterceptor.USER_EXCEPTION;

import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleGroup;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;

public class EllipsoidCalcController {

    @FXML
    private TextField phiStart;
    
    @FXML
    private ImageView diagramImage;
    
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
    	Image image = new Image(EllipsoidCalcMain.class.getResourceAsStream("/resources/PcoordinatesImage.jpg"));
		diagramImage.setImage(image);
    }
    
    public void convertDecimal(ActionEvent e){
    	TextField temp = (TextField) e.getSource();
    	try {
			temp.setText(String.format("%f", EllipsoidCalcUtility.convertToDecimal(temp.getText())));
		} catch (InvalidUserInputException e1) {
			warningLabel.setText(String.format("Non-number '%s' detected! Fix this to proceed", e1.getMessage()));
		}
    	
    }
    
    public void calculateButton(){
    	if(ensureAllEntriesLogged()){
    		return;
    	}
    	
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
			warningLabel.setText(String.format("Non-number '%s' detected", e1.getMessage()));
		}
    	try {
			double thetaBegin = EllipsoidCalcUtility.convertThetaToRadians(startAngle, inDegrees);
			double thetaEnd = EllipsoidCalcUtility.convertThetaToRadians(endAngle, inDegrees);
			double phiAngleEnd = EllipsoidCalcUtility.convertThetaToRadians(zAxisAngleEnd, inDegrees);
			double phiAngleStart = EllipsoidCalcUtility.convertThetaToRadians(zAxisAngleStart, inDegrees);
			volumeAnswer.setText(String.format("%f", EllipsoidCalcUtility.getVolume(thetaBegin, thetaEnd, phiAngleEnd, phiAngleStart,
					aAxis, bAxis, cAxis)));
		} catch (InvalidUserInputException e) {
			warningLabel.setText(String.format("Non number '%s' Detected", e.getMessage()));
		}
    	
    }
    
    private boolean ensureAllEntriesLogged(){
    	boolean incompleteInput = false;
    	
    	if(phiEnd.getText().equals("")){
    		warningLabel.setText("All fields must contain a number");
    		incompleteInput = true;
    	}
    	
    	if(phiStart.getText().equals("")){
    		warningLabel.setText("All fields must contain a number");
    		incompleteInput = true;
    	}
    	
    	if(thetaStart.getText().equals("")){
    		warningLabel.setText("All fields must contain a number");
    		incompleteInput = true;
    	}
    	
    	if(thetaEnd.getText().equals("")){
    		warningLabel.setText("All fields must contain a number");
    		incompleteInput = true;
    	}
    	
    	if(a.getText().equals("")){
    		warningLabel.setText("All fields must contain a number");
    		incompleteInput = true;
    	}
    	
    	if(b.getText().equals("")){
    		warningLabel.setText("All fields must contain a number");
    		incompleteInput = true;
    	}
    	
    	if(c.getText().equals("")){
    		warningLabel.setText("All fields must contain a number");
    		incompleteInput = true;
    	}
    	
    	return incompleteInput;
    }
    
    
    
}