package analyzer.databaseRunners;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.text.DecimalFormat;

import analyzer.Utilities.Utilities;
import analyzer.variantInfo.Variant;

public class pdbBlastRunner {

    String tempfaaPath;
    String tempoutPath;

    
    public pdbBlastRunner(String outputPath) throws IOException{
        this.tempfaaPath = outputPath + "temp.faa";
        this.tempoutPath = outputPath + "temp.out";

        }
    
    public void runPDBBlast(Variant var) throws IOException {
    	
    	//System.out.println(var.toString());
    	
    	if (var.getCDSList().size() == 1) {
    		buildPDBQuery(var, 0); 	
        	try {
    			runPDBBlastCommand();
    		} catch (Exception e) {
    			e.printStackTrace();
    		}
        	extractPDBBlastResults(var, 0);
    	} else if (var.getCDSList().size() != 0){
    		for (int i = 0; i < var.getCDSList().size(); ++i) {
    			buildPDBQuery(var, i); 	
            	try {
        			runPDBBlastCommand();
        		} catch (Exception e) {
        			e.printStackTrace();
        		}
            	extractPDBBlastResults(var, i);
    		}
    	}
    	
    	
    }

    public void buildPDBQuery(Variant var, Integer num) throws IOException {
    
    	FileWriter tempFAA = new FileWriter(tempfaaPath);
    	StringBuilder fasta = new StringBuilder('>' + var.getChr() + ':' + var.getPos() + '\n' +
    											var.getCDSList().get(num).getOriginalProtein().substring(0, var.getCDSList().get(num).getOriginalProtein().length() - 1) //take off *
    											+ '\n');
    	tempFAA.write(fasta.toString());
    	tempFAA.close();
    }
    
    public void runPDBBlastCommand() throws Exception {
    	System.out.println(Utilities.GREEN+"Running blastp on PDB"+ Utilities.RESET);
    	try {
    		String[] call = new String[]{"blastp", "-query", tempfaaPath, "-db", "pdb",
    									 "-out", tempoutPath, "-evalue", ".005", "-outfmt", "6 sseqid qstart qend length evalue stitle"};
    		
    		ProcessBuilder pb = new ProcessBuilder(call);

    		Process p = pb.start();
    		p.waitFor();
    		//System.out.println(Utilities.getProcessOutput(p));
    		String error = Utilities.getProcessError(p);
    		p.destroy();
    		
    		Files.deleteIfExists(new File(tempfaaPath).toPath());
    		
    		if(!error.isEmpty()){
                throw new Exception("blastp (with pdb db) threw the following error: " + error);
            }
    	} catch (IOException e) {
    		e.printStackTrace();
    	} catch (InterruptedException e) {
    		e.printStackTrace();
    	}
    }
    
    public void extractPDBBlastResults(Variant var, Integer num) throws IOException {
    	String thisLine = null;
    	
    	BufferedReader pdbResults = new BufferedReader(new FileReader(new File(tempoutPath)));
    	while ((thisLine = pdbResults.readLine()) != null) {
    		
    		
    		//4j4j_A	8	180	180	4e-28	4j4j_A mol:protein length:209  DNA dC->dU-editing enzyme APOBEC-3F

    		String[] splitLine = thisLine.split("\t");
    		
    		Integer cddStart = Integer.parseInt(splitLine[1]);
    		Integer cddEnd = Integer.parseInt(splitLine[2]);
    		Double percentLost;
            Double position = Math.ceil(var.WithinGenePos.get(num) / 3.0); // convert to the protein position
            
    		DecimalFormat df = new DecimalFormat("#.##");
			df.setRoundingMode(RoundingMode.CEILING);
			
    		if (cddStart >= position) { // starts after the lost splice site
    			percentLost = 100.0;

    		} else if (cddEnd >= position) { // variant within this domain
    			Double totalDomainLength = (double) (cddEnd - cddStart + 1);
    			Double lostAmount = (double) (cddEnd - position + 1);
    			percentLost = lostAmount / totalDomainLength * 100;
    			

    		} else {
    			percentLost = 0.0;
    		}

    		
    		if (!percentLost.equals(0.0)) {
    			
    			StringBuilder outCD = new StringBuilder();
        		outCD.append(var.getCDSList().get(num).getTransName() + ':' + var.getCDSList().get(num).cDot + "\t" + 
        					 splitLine[0] + '\t' + // gnl|CDD|306940
        					 df.format(percentLost) + "%\t" + // percentLost
        					 splitLine[4] + '\t' + // e-val for match
        					 splitLine[5] + '\n'
        				);
        	
    		var.PDBList.add(outCD.toString());
        	
       		}
    		
         }  
    	
    	
		Files.deleteIfExists(new File(tempoutPath).toPath());
    	
    	pdbResults.close();
    }
    
}