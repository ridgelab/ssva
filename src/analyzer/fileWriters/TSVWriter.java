package analyzer.fileWriters;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.lang.StringBuilder;

import analyzer.transcriptInfo.CDS;
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
    	
    	
    }
    
    public void writeVariant(Variant var) throws IOException {
    	
    	/* private String Chr; // chr21
    private Integer Pos; // 1827643
    private String Ref; // A
    private String Alt; // T
    private String homhet; // hom het
    private String spliceInfo; // NONE(dist=NONE),MIR3648(dist=414230)
    private ArrayList<String> Annotations;
    private ArrayList<CDS> transcripts;
    private String GeneName;
    private ArrayList<Double> OriginalMesScores;
    private ArrayList<Double> VariantMesScores;
    private ArrayList<Double> percentDiffList;*/
    	
    	StringBuilder variantTSV = new StringBuilder();
    	variantTSV.append(var.getChr() + '\t' +
    					  var.getPos() + '\t' +
    					  var.getRef() + '\t' +
    					  var.getAlt());
    	
    	file.write(variantTSV.toString());
    }
	
    public String getPath(){
        return this.path;
    }

    public void close() throws IOException {
        this.file.close();
    }
    
}
