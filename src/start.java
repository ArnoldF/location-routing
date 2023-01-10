import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Random;

import org.json.simple.JSONArray;
import org.json.simple.JSONObject;

import MDVRP.Stopwatch;

public class start {
	
	public static void main(String[] args) {
		
		String path = readInputArgs(args);
		String outputFile = "output";//args[1];
		
		LRP lrp = new LRP(1);
		
		lrp.loadInstance(path);	
		Stopwatch stopwatch = new Stopwatch();
		stopwatch.start();
		lrp.solve(new Integer[]{100,10,3,1});
		stopwatch.stop();
		//System.out.println(stopwatch.elapsedTimeMillis() / 1000.0);
		//System.out.println();
		writeToFile(lrp, stopwatch.elapsedTimeMillis(), outputFile);
		//lrp.writeToFile();	
	}

	private static String readInputArgs(String[] args) {
		boolean hasMinCapacity = false;
		for (String input: args) {
			if (input.substring(input.lastIndexOf('.') + 1).equals("json")) {
				return input;
			}
		}
		return "";
	}
	
	/*private static void readInputArgs(String[] args, LRPScenarios lrp) {
		boolean hasMinCapacity = false;
		for (String input: args) {
			if (input.charAt(0) == '-') {
				if (input.contains("v")) {
					lrp.setValidate();
				}
				if (input.contains("m")) {
					hasMinCapacity = true;
				}
			}
			else if (lrp.isValidation() && input.substring(input.lastIndexOf('.') + 1).equals("sol")) {
				lrp.readSolution(input);
			}
			else if (hasMinCapacity && input.matches("-?\\d+")) {
				lrp.setMinCapacity(Integer.parseInt(input));
			}
			else if (input.substring(input.lastIndexOf('.') + 1).equals("json")) {
				lrp.addScenario(input);	
			}
		}
	}*/
	
	private static ArrayList<String> listFilesForFolder(final File folder) 
	{
		ArrayList<String> files = new ArrayList<String>();
	    for (final File fileEntry : folder.listFiles()) {
	        if (fileEntry.isDirectory()) {
	            files.addAll(listFilesForFolder(fileEntry));
	        } else {
	        	files.add(fileEntry.getName());
	        }
	    }
	    return files;
	}
	
	
	public static void writeToFile(LRP lrp, long time, String outputFile) {
        String path = "." + File.separator +"output" + File.separator + outputFile+ ".sol";

	    try (FileWriter file = new FileWriter(path, true)) {
	        file.write(lrp.remainingConfigurations.first().costs + " " + time/1000 + " " + lrp.remainingConfigurations.first().openDepots.size() + " " + lrp.upperBoundOpenDepots + "\r\n");

	    } catch (IOException e) {
	    }
	}


	
}
