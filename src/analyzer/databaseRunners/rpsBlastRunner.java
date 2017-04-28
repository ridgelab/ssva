package analyzer.databaseRunners;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;

import analyzer.Utilities.Utilities;
import analyzer.variantInfo.Variant;

public class rpsBlastRunner {

    String tempfaaPath;
    String tempoutPath;

    
    public rpsBlastRunner(String outputPath) throws IOException{
        this.tempfaaPath = outputPath + "temp.faa";
        this.tempoutPath = outputPath + "temp.out";

        }
    
    public void runRPSBlast(Variant var) throws IOException {
    	
    	buildRPSQuery(var); 	
    	runRPSBlastCommand();
    	extractRPSBlastResults();
    	
    }

    public void buildRPSQuery(Variant var) throws IOException {
    
    	FileWriter tempFAA = new FileWriter(tempfaaPath);
    	StringBuilder fasta = new StringBuilder('>' + var.getChr() + ':' + var.getPos() + '\n' +
    											var.getCDSList().get(0).getOriginalProtein() + '\n');
    	tempFAA.write(fasta.toString());
    	tempFAA.close();
    }
    
    public void runRPSBlastCommand() {
    	System.out.println(Utilities.GREEN+"Running rpsblast to find Conserved Domains"+ Utilities.RESET);
    	try {
    		String[] call = new String[]{"rpsblast", "-query", tempfaaPath, "-db Cdd",
    									 "-out", tempoutPath, "-evalue .05", "-outfmt '6 OPTIONS'"};
    		ProcessBuilder pb = new ProcessBuilder(call);

    		Process p = pb.start();
    		p.waitFor();
    		System.out.println(Utilities.getProcessOutput(p));
    		System.out.println(Utilities.getProcessError(p));
    		p.destroy();
    		
    		Files.deleteIfExists(new File(tempfaaPath).toPath());
    		
    	} catch (IOException e) {
    		e.printStackTrace();
    	} catch (InterruptedException e) {
    		e.printStackTrace();
    	}
    }
    
    public void extractRPSBlastResults() throws IOException {
    	String thisLine = null;
    	
    	BufferedReader rpsResults = new BufferedReader(new FileReader(new File(tempoutPath)));
    	while ((thisLine = rpsResults.readLine()) != null) {
            System.out.println(thisLine);
         }  
    	
		Files.deleteIfExists(new File(tempoutPath).toPath());
    	
    	rpsResults.close();
    	
    }
    
}