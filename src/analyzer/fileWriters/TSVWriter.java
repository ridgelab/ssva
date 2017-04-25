package analyzer.fileWriters;

import java.io.FileWriter;
import java.io.IOException;
import java.lang.StringBuilder;

import analyzer.Utilities.Utilities;
import analyzer.variantInfo.Variant;

/**
 * Created by CS on 4/25/17.
 */
public class TSVWriter {

    private String path;
    FileWriter file;

    public TSVWriter(String tsvPath) throws IOException{
    	path = tsvPath;
    	file = new FileWriter(tsvPath);
    	
    	file.write("CHR\tPOS\tREF\tALT\tGENE\tGERP2\t1000GEN\tEXAC\tWT_MESSCORE\tVAR_MESSCORE\t%DIFF\tSPLICEINFO\n");
    }
    
    public void writeVariant(Variant var) throws IOException {
    	
    	/* private String Chr; // chr21
    private Integer Pos; // 1827643
    private String Ref; // A
    private String Alt; // T
    private String homhet; // hom het
    private String spliceInfo; // NONE(dist=NONE),MIR3648(dist=414230)
    public ArrayList<String> Annotations;
    private ArrayList<CDS> transcripts;
    private String GeneName;
    private ArrayList<Double> OriginalMesScores;
    private ArrayList<Double> VariantMesScores;
    private ArrayList<Double> percentDiffList;*/
    	
    	StringBuilder variantTSV = new StringBuilder();
    	variantTSV.append(var.getChr() + '\t' +
    					  var.getPos() + '\t' +
    					  var.getRef() + '\t' +
    					  var.getAlt() + '\t' +
    					  var.getGeneName() + '\t'+
    					  var.Annotations.get(0) + '\t' +
    					  var.Annotations.get(1) + '\t' +
    					  var.Annotations.get(2) + '\t' +
    					  var.getOriginalMesScores().get(0) + '\t');
    	
		for	(double MESScore : var.getVariantMesScores()) {
    		variantTSV.append(Double.toString(MESScore) +';');
            System.out.println(Utilities.GREEN + "MESScore: " + Utilities.RESET + Double.toString(MESScore));

    	}
		
		variantTSV.append('\t');
		
		for	(double DiffScore : var.getPercentDiffList()) {
    		variantTSV.append(Double.toString(DiffScore) +';');
            System.out.println(Utilities.GREEN + "DiffScore: " + Utilities.RESET + Double.toString(DiffScore));

    	}
		
		variantTSV.append('\t' + var.getSpliceInfo() + '\n');
    	
    	file.write(variantTSV.toString());
    }
	
    public String getPath(){
        return this.path;
    }

    public void close() throws IOException {
        this.file.close();
    }
    
}
