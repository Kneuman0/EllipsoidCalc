package biz.personalAcademics.ellipsoidCalc;

import java.text.DecimalFormat;

import biz.personalAcademics.ellipsoid.Ellipsoid;
import biz.personalAcademics.ellipsoid.customExceptions.InvalidUserInputException;
import javafx.application.Platform;
import javafx.event.ActionEvent;
import javafx.fxml.FXML;
import javafx.scene.control.Alert;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.TextField;
import javafx.scene.control.ToggleGroup;
import javafx.scene.image.Image;
import javafx.scene.image.ImageView;
import javafx.scene.layout.AnchorPane;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundImage;
import javafx.scene.layout.BackgroundPosition;
import javafx.scene.layout.BackgroundRepeat;
import javafx.scene.layout.BackgroundSize;

public class EllipsoidCalcController {

	@FXML
	private TextField 
			sampleSizeTextBox,
			phiStart,
			phiEnd,
			thetaStart,
			thetaEnd,
			a,
			b,
			c;

	@FXML
	private Label 
			warningLabel,
			volumeAnswer;

	@FXML
	private RadioButton 
			radiansRadio,
			degreesRadio;

	@FXML
	private ImageView diagramImage;

	@FXML
	private ToggleGroup angle;

	@FXML
	private AnchorPane anchorPane;

	public void initialize() {
		Image image = new Image(
				EllipsoidCalcMain.class
						.getResourceAsStream("/resources/PcoordinatesImage.jpg"));
		diagramImage.setImage(image);
		setBackground();
		presentProgramInformationToUser();
	}

	public void convertDecimal(ActionEvent event) {
		TextField temp = (TextField) event.getSource();
		try {
			temp.setText(String.format(
					"%f",
					Ellipsoid.convertToDecimal(temp.getText().replaceAll(" ",
							""))));
		} catch (InvalidUserInputException e) {
			warningLabel.setText(String.format(
					"Non-number '%s' detected! Fix this to proceed",
					e.getMessage()));
		}

	}

	public void calculateButton() {
		warningLabel.setText("");
		volumeAnswer.setText("");
		if (ensureAllEntriesLogged()) {
			return;
		}

		String startAngle = thetaStart.getText().replaceAll(" ", "");
		String endAngle = thetaEnd.getText().replaceAll(" ", "");
		String zAxisAngleEnd = phiEnd.getText().replaceAll(" ", "");
		String zAxisAngleStart = phiStart.getText().replaceAll(" ", "");
		boolean inDegrees = degreesRadio.isSelected();
		double aAxis = 0;
		double bAxis = 0;
		double cAxis = 0;
		try {
			aAxis = Double.parseDouble(a.getText().replaceAll(" ", ""));
			bAxis = Double.parseDouble(b.getText().replaceAll(" ", ""));
			cAxis = Double.parseDouble(c.getText().replaceAll(" ", ""));
		} catch (NumberFormatException e1) {
			warningLabel.setText(String.format("Non-number '%s' detected",
					e1.getMessage()));
			return;
		}

		Ellipsoid ellip = null;
		Thread calculation = null;

		try {
			double thetaBegin = Ellipsoid.convertThetaToRadians(startAngle,
					inDegrees);
			double thetaEnd = Ellipsoid.convertThetaToRadians(endAngle,
					inDegrees);
			double phiAngleEnd = Ellipsoid.convertThetaToRadians(zAxisAngleEnd,
					inDegrees);
			double phiAngleStart = Ellipsoid.convertThetaToRadians(
					zAxisAngleStart, inDegrees);

			ellip = new Ellipsoid(thetaBegin, thetaEnd, phiAngleStart,
					phiAngleEnd, aAxis, bAxis, cAxis);
			
			/*
			 * Creates a new thread to run calculation. Inner class accepts
			 * an Ellipsoid object in its constructor
			 */
			calculation = new Thread(new ExecuteCalculation(ellip));

		} catch (InvalidUserInputException e) {
			warningLabel.setText(String.format("Non number '%s' Detected",
					e.getMessage()));
			volumeAnswer.setText("");
			// terminate is invalid user input is detected
			return;
		}
		
		/*
		 * Starts the calculation in a new thread. This informs the user
		 * that the calculation is in progress and then displays the result
		 */
		calculation.start();

	

	}
	

	private boolean ensureAllEntriesLogged() {
		boolean incompleteInput = false;

		if (phiEnd.getText().equals("")) {
			warningLabel.setText("All fields must contain a number");
			incompleteInput = true;
		}

		if (phiStart.getText().equals("")) {
			warningLabel.setText("All fields must contain a number");
			incompleteInput = true;
		}

		if (thetaStart.getText().equals("")) {
			warningLabel.setText("All fields must contain a number");
			incompleteInput = true;
		}

		if (thetaEnd.getText().equals("")) {
			warningLabel.setText("All fields must contain a number");
			incompleteInput = true;
		}

		if (a.getText().equals("")) {
			warningLabel.setText("All fields must contain a number");
			incompleteInput = true;
		}

		if (b.getText().equals("")) {
			warningLabel.setText("All fields must contain a number");
			incompleteInput = true;
		}

		if (c.getText().equals("")) {
			warningLabel.setText("All fields must contain a number");
			incompleteInput = true;
		}

		return incompleteInput;
	}

	private void setBackground() {
		Image logo = new Image(
				EllipsoidCalcMain.class
						.getResourceAsStream("/resources/mathBackground.jpg"));
		BackgroundSize logoSize = new BackgroundSize(600, 400, false, false,
				true, true);
		BackgroundImage image = new BackgroundImage(logo,
				BackgroundRepeat.NO_REPEAT, BackgroundRepeat.NO_REPEAT,
				BackgroundPosition.CENTER, logoSize);
		Background background = new Background(image);
		anchorPane.setBackground(background);
	}

	private void presentProgramInformationToUser() {
		String sphericalCoordinatesWarning = "This program requires that you have a basic understanding\nof"
				+ " ellipsoids and spherical coordinates";
		String axisWarning = "the a axis is associated with the x axis,\nthe b axis is associated with the y axis,\n"
				+ "and the c axis is associated with the z axis.";
		Alert alert = new Alert(Alert.AlertType.INFORMATION);
		alert.setTitle("Usage Information");
		alert.setHeaderText(null);
		alert.setContentText(String.format("%s\n\n%s",
				sphericalCoordinatesWarning, axisWarning));
		alert.showAndWait();
	}
	
	private class ExecuteCalculation implements Runnable{
		
		private Ellipsoid ellip;
		
		public ExecuteCalculation(Ellipsoid ellip) {
			this.ellip = ellip;
		}
		
		@Override
		public void run() {
			
			try {
				if (sampleSizeTextBox.getText().equals("")) {
					Platform.runLater(new UpdateLabel(warningLabel, "Calculating... Please wait"));
					Platform.runLater(new UpdateLabel(volumeAnswer, ellip.toString()));
					Platform.runLater(new UpdateLabel(warningLabel, ""));
				} else {
					String sampleString = sampleSizeTextBox.getText().replaceAll("[,_]", "");
					int sampleSize = Integer.parseInt(sampleString);
					
					Platform.runLater(new UpdateLabel(warningLabel, "Calculating... Please wait"));
					Platform.runLater(new UpdateLabel(volumeAnswer, ellip.toString(sampleSize)));
					Platform.runLater(new UpdateLabel(warningLabel, ""));
				}
			} catch (NumberFormatException | InvalidUserInputException e) {
				DecimalFormat sampleForm = new DecimalFormat("#,###,###");
				Platform.runLater(new UpdateLabel(warningLabel, String
						.format("%s is not an acceptable sample size, sample must be at least %s",
								e.getMessage(), sampleForm.format(Ellipsoid.MIN_SAMPLE_SIZE))));
				return;
			}
			
			
		}
		
	}
	
	private class UpdateLabel implements Runnable{
		
		private Label label;
		private String message;
		
		public UpdateLabel(Label label, String message) {
			this.label = label;
			this.message = message;
		}
		
		@Override
		public void run() {
			label.setText(message);
		}
		
	}

}