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
    FileWriter file;

    public TSVWriter(String tsvPath) throws IOException{
    	path = tsvPath;
    	file = new FileWriter(tsvPath);
    	
    	file.write("#CHR\tPOS\tREF\tALT\tGENE\tGERP2\t1000GEN\tEXAC\tWT_MESSCORE\tVAR_MESSCORE\t%DIFF\tSPLICEINFO\n");
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

    	/* ROUNDING
    	 * DecimalFormat df = new DecimalFormat("#.####");
			df.setRoundingMode(RoundingMode.CEILING);
			for (Number n : Arrays.asList(12, 123.12345, 0.23, 0.1, 2341234.212431324)) {
    		Double d = n.doubleValue();
    		System.out.println(df.format(d));
			}
    	 */
    	
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
    					  var.Annotations.get(1) + '\t' +
    					  var.Annotations.get(2) + '\t'
    					  );
    	
    	
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
    		variantTSV.append(df.format(var.getPercentDiffList().get(0).doubleValue()));
    	} else if (var.getPercentDiffList().size() == 0){
    		variantTSV.append("NA\t");
    	} else {
    		for	(double DiffScore : var.getPercentDiffList()) {
    			variantTSV.append(df.format(DiffScore) +';');
    		}
    	}
		
		variantTSV.append('\t' + var.getSpliceInfo() + '\n');
    	
    	file.write(variantTSV.toString());
    }
    
    public void writeConservedDomains(Variant var) throws IOException {
    	file.write("#CDDid\tSTART\t%LOST\tE-VAL\tINFO\n");
    	
    	StringBuilder variantTSV = new StringBuilder();
    	variantTSV.append(var.getChr() + '\t' +
    					  var.getPos() + '\t' +
    					  var.getRef() + '\t' +
    					  var.getAlt() + '\t' +
    					  var.getGeneName() + '\t'+
    					  var.Annotations.get(0) + '\t' +
    					  var.Annotations.get(1) + '\t' +
    					  var.Annotations.get(2));    	

    	
    	file.write(variantTSV.toString());
    }
	
    public String getPath(){
        return this.path;
    }

    public void close() throws IOException {
        this.file.close();
    }
    
}
