package analyzer.fileWriters;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.StringBuilder;
import java.math.RoundingMode;
import java.text.DecimalFormat;

//import analyzer.Utilities.Utilities;
import analyzer.variantInfo.Variant;

/**
 * Created by CS on 4/25/17.
 */
public class TSVWriter {

    private String path;
    public FileWriter file;
    public String build;

    public TSVWriter(String tsvPath, String build) throws IOException{
    	path = tsvPath;
    	file = new FileWriter(tsvPath);

    	this.build = build;
    	if (build.equals("hg19")) {
        	file.write("#CHR\tPOS\tREF\tALT\tGENE\t1000GEN\tEXAC\tGERP2\tWT_MESSCORE\tVAR_MESSCORE\t%DIFF\tSPLICEINFO\n");
    	} else if (build.equals("hg38")) {
        	file.write("#CHR\tPOS\tREF\tALT\tGENE\t1000GEN\tEXAC\tWT_MESSCORE\tVAR_MESSCORE\t%DIFF\tSPLICEINFO\n");
    	}
    }
    
    public void writeVariant(Variant var) throws IOException {

    	// rounding
    	DecimalFormat df = new DecimalFormat("#.####");
    	df.setRoundingMode(RoundingMode.CEILING);
    	
    	
    	StringBuilder variantTSV = new StringBuilder();
    	
    	variantTSV.append(var.getChr() + '\t' +
    					  var.getPos() + '\t' +
    					  var.getRef() + '\t' +
    					  var.getAlt() + '\t' +
    					  var.getGeneName() + '\t'+
    					  var.Annotations.get(0) + '\t' +
    					  var.Annotations.get(1) + '\t');
    	
    	if (build.equals("hg19")) {
    		variantTSV.append( var.Annotations.get(2) + '\t');
    	}
    	
    	if (var.getOriginalMesScores().size() != 0) {
    		variantTSV.append(df.format(var.getOriginalMesScores().get(0)) + '\t');
    	} else {
    		variantTSV.append("NA\t");
    	}
    	
    	if (var.getVariantMesScores().size() == 1) {
    		variantTSV.append(df.format(var.getVariantMesScores().get(0).doubleValue()));
    	} else if (var.getVariantMesScores().size() == 0){
    		variantTSV.append("NA\t");
    	} else {
    		for	(Double MESScore : var.getVariantMesScores()) {
        		variantTSV.append(df.format(MESScore) +';');
        	}
    	}
		
		variantTSV.append('\t');
		
		if (var.getPercentDiffList().size() == 1) {
    		variantTSV.append(df.format(var.getPercentDiffList().get(0).doubleValue()) + '%');
    	} else if (var.getPercentDiffList().size() == 0){
    		variantTSV.append("NA\t");
    	} else {
    		for	(double DiffScore : var.getPercentDiffList()) {
    			variantTSV.append(df.format(DiffScore) +';');
    		}
    	}
		
		variantTSV.append('\t' + var.getSpliceInfo() + '\n');
    	
    	file.write(variantTSV.toString());
    	if (var.ConservedDomains.size() != 0) {
        	this.writeConservedDomains(var);
    	}
    	/*if (var.PDBList.size() != 0) {
        	this.writePDB(var);
    	}*/

    }
    
    public void writeConservedDomains(Variant var) throws IOException {
    	file.write("\t#TRANSCRIPT\tCDDid\t%LOST\tE-VAL\tINFO\n");

    	for (String out : var.ConservedDomains) {
        	file.write('\t' + out);
    	}
    }
    
    /*public void writePDB(Variant var) throws IOException {
    	file.write("\t#TRANSCRIPT\tPDBid\t%LOST\tE-VAL\tINFO\n");

    	for (String out : var.PDBList) {
        	file.write('\t' + out);
    	}
    }*/
	
    public String getPath(){
        return this.path;
    }

    public void close() throws IOException {
        this.file.close();
    }
   
}
