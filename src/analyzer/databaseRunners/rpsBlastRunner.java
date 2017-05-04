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
    	try {
			runRPSBlastCommand();
		} catch (Exception e) {
			e.printStackTrace();
		}
    	extractRPSBlastResults(var);
    	
    }

    public void buildRPSQuery(Variant var) throws IOException {
    
    	FileWriter tempFAA = new FileWriter(tempfaaPath);
    	StringBuilder fasta = new StringBuilder('>' + var.getChr() + ':' + var.getPos() + '\n' +
    											var.getCDSList().get(0).getOriginalProtein().substring(0, var.getCDSList().get(0).getOriginalProtein().length() - 1) //take off *
    											+ '\n');
    	tempFAA.write(fasta.toString());
    	tempFAA.close();
    }
    
    public void runRPSBlastCommand() throws Exception {
    	System.out.println(Utilities.GREEN+"Running rpsblast to find Conserved Domains"+ Utilities.RESET);
    	try {
    		String[] call = new String[]{"rpsblast", "-query", tempfaaPath, "-db", "Cdd",
    									 "-out", tempoutPath, "-evalue", ".01", "-outfmt", "6 sseqid qstart qend length evalue stitle"};
    		
    		ProcessBuilder pb = new ProcessBuilder(call);

    		Process p = pb.start();
    		p.waitFor();
    		//System.out.println(Utilities.getProcessOutput(p));
    		String error = Utilities.getProcessError(p);
    		p.destroy();
    		
    		Files.deleteIfExists(new File(tempfaaPath).toPath());
    		
    		if(!error.isEmpty()){
                throw new Exception("rpsblast threw the following error: " + error);
            }
    	} catch (IOException e) {
    		e.printStackTrace();
    	} catch (InterruptedException e) {
    		e.printStackTrace();
    	}
    }
    
    public void extractRPSBlastResults(Variant var) throws IOException {
    	String thisLine = null;
    	
    	BufferedReader rpsResults = new BufferedReader(new FileReader(new File(tempoutPath)));
    	while ((thisLine = rpsResults.readLine()) != null) {
            System.out.println(thisLine);
         }  
    	
		Files.deleteIfExists(new File(tempoutPath).toPath());
    	
    	rpsResults.close();
    	
    }
    
}