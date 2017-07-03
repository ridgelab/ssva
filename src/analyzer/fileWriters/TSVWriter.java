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
    public String rps;


    public TSVWriter(String tsvPath, String build, String rps) throws IOException{
    	path = tsvPath;
    	file = new FileWriter(tsvPath);
    	this.rps = rps;

    	this.build = build;
    	if (build.equals("hg19")) {
    		if (rps.equals("true")) {
            	file.write("#CHR\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\t1000GEN\tEXAC\tdbscSNV\tGERP2\tWT_MESSCORE\tVAR_MESSCORE\t%DIFF\tSPLICEINFO\tCDD(TRANSCRIPT//CDDid//%LOST//E-VAL//INFO)\n");
    		} else {
            	file.write("#CHR\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\t1000GEN\tEXAC\tdbscSNV\tGERP2\tWT_MESSCORE\tVAR_MESSCORE\t%DIFF\tSPLICEINFO\n");
    		}
    	} else if (build.equals("hg38")) {
    		if (rps.equals("true")) {
            	file.write("#CHR\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\t1000GEN\tEXAC\tdbscSNV\tWT_MESSCORE\tVAR_MESSCORE\t%DIFF\tSPLICEINFO\tCDD(TRANSCRIPT//CDDid//%LOST//E-VAL//INFO)\n");
    		} else {
            	file.write("#CHR\tPOS\tREF\tALT\tGENE\tTRANSCRIPT\t1000GEN\tEXAC\tdbscSNV\tWT_MESSCORE\tVAR_MESSCORE\t%DIFF\tSPLICEINFO\n");
    		}
    	}
    }
    
    public void flushIt() throws IOException {
    	file.flush();
    }
    
    public void writeVariant(Variant var) throws IOException {

    	for (int i = 0; i < var.getCDSList().size(); ++i) {
    		// rounding
        	DecimalFormat df = new DecimalFormat("#.####");
        	df.setRoundingMode(RoundingMode.CEILING);	
        	
        	StringBuilder variantTSV = new StringBuilder();
        	
        	variantTSV.append(var.getChr() + '\t' +
        					  var.getPos() + '\t' +
        					  var.getRef() + '\t' +
        					  var.getAlt() + '\t' +
        					  var.getGeneName() + '\t'+
        					  var.getCDSList().get(i).transName + '\t' +
        					  var.Annotations.get(0) + '\t' +
        					  var.Annotations.get(1) + '\t' +
        					  var.Annotations.get(2) + '\t');
        	
        	if (build.equals("hg19")) {
        		variantTSV.append( var.Annotations.get(3) + '\t');
        	}
        	
        	if (var.getOriginalMesScores().size() != 0) {
        		variantTSV.append(df.format(var.getOriginalMesScores().get(0)) + '\t');
        	} else {
        		variantTSV.append("NA\t");
        	}
        	        	
        	if (var.getVariantMesScores().size() != 0) {
        		variantTSV.append(df.format(var.getVariantMesScores().get(i).doubleValue()) + '\t');
        	} else {
        		variantTSV.append("NA\t");
        	}
    		    		
    		if (var.getPercentDiffList().size() != 0) {
        		variantTSV.append(df.format(var.getPercentDiffList().get(i).doubleValue()) + '%');
        	} else {
        		variantTSV.append("NA\t");
        	} 
    		
    		variantTSV.append('\t' + var.getCDSList().get(i).getcDot() + '\t');

			file.write(variantTSV.toString());


    		if (rps.equals("true")) {
            	if (var.ConservedDomains.size() != 0) {
                	this.writeConservedDomains(var, var.getCDSList().get(i).transName);
            	}
    		} else {
    	    	file.write('\n');
    		}
        	
        	
        	
    	}
    	/*if (var.PDBList.size() != 0) {
        	this.writePDB(var);
    	}*/

    }
    
    public void writeConservedDomains(Variant var, String transName) throws IOException {
    	//file.write("\t#TRANSCRIPT\tCDDid\t%LOST\tE-VAL\tINFO\n");

    	boolean first = true;
    	
    	for (String out : var.ConservedDomains) {
    		if (out.contains(transName)) {
    			if (first) {
        			file.write(out);
    			} else {
    				file.write(out);
    			}
    		}
    	}
    	file.write('\n');
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
