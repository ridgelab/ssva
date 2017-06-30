package analyzer;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import analyzer.RefSeq.PullRegionsFromRef;
import analyzer.RefSeq.RefSeqParser;
import analyzer.Utilities.Utilities;
import analyzer.annovarParsers.GeneParser;
import analyzer.annovarParsers.GeneralAnnotationParser;
import analyzer.databaseRunners.AnnovarRunner;
import analyzer.databaseRunners.MESRunner;
//import analyzer.databaseRunners.pdbBlastRunner;
import analyzer.databaseRunners.rpsBlastRunner;
import analyzer.fileWriters.TSVWriter;
//import analyzer.fileWriters.VCFWriter;
import analyzer.fileWriters.annovarWriter;
import analyzer.variantInfo.Variant;
//import htsjdk.variant.variantcontext.VariantContextBuilder;
//import analyzer.transcriptInfo.CDS;
import net.sourceforge.argparse4j.inf.Namespace;

/*
 * This class is used to run the splice detection algorithm. 
 */
public class SpliceRunner {

    private String annovar;                // Path to annovar 
    private String input;                  // Input file in vcf file 
    private String human;                  // Path to the humandb to be used by annovar 
    private String outputFolder;           // Directory to save outputfiles to
    private String ref;                    // Path to the directory containing the UCSC reference genome
    private String refSeq;                 // Path to the file that contains the refseq data 
    private TreeMap<String,Variant> vars;  // A map of the splicing variants. keys = chromsome:position, values = splicing variant. 
    private String SamtoolsPath;           // Path to Samtools 
    private String MaxEntPath;             // Path to the algorithm directory 
    private String build;				   // Build version of genome (default: hg19)
    private String eval;				   // rpsblast e-value cut-off
    private String rpsblast;			   // whether to to rpsblast or not



    //-------------------------------------------------------------------------------------
    // Constructor method 
    //-------------------------------------------------------------------------------------
    public SpliceRunner(Namespace res){

        this.annovar = res.getString("Annovar");
        this.input = res.getString("Input");
        this.human = res.getString("human");
        this.outputFolder = res.getString("Output");
        if (!this.outputFolder.endsWith("/")){
            StringBuilder sb = new StringBuilder(this.outputFolder+"/");
            this.outputFolder = sb.toString();
        }
        this.ref = res.getString("Ref");
        this.refSeq = res.getString("analyzer/RefSeq");
        this.build = res.getString("build");
        if (this.build.equals("hg38") && this.refSeq.equals("RefSeq_hg19.txt")){
        	this.refSeq = "RefSeq_hg38.txt";
        }
        this.SamtoolsPath = res.getString("Samtools");
        this.MaxEntPath = res.getString("MaxEntPath");
        this.eval = res.getString("eval");
        this.rpsblast = res.getString("rpsblast");

    }

    //-------------------------------------------------------------------------------------
    // Method used to run the splice site algorithm
    //-------------------------------------------------------------------------------------
    public void run() throws Exception{
        String varfunc = convertAnnovar(); // Convert input VCF to annovar format and run annovar on it. (Returns the file name of the file, containing annotated variants, created by annovar)
        String newFile = geneBasedAnnotation(varfunc); // SETS VARS / Parse the annotated file. (Returns the file name of a "splice file" containing splice variants) 
        runAnnotations(newFile); // Runs annovar 3 times to get 3 different conserved scores. New information is added to the "vars" map.
        RefSeqParser rsp = new RefSeqParser(this.refSeq); //Parses the refseq data file. (2 files created) 
        PullRegionsFromRef prfr = new PullRegionsFromRef(ref,SamtoolsPath);  //hg19 / Sam tools
        Iterator<Map.Entry<String,Variant>> iter = this.vars.entrySet().iterator();
        rpsBlastRunner rpsRunner = new rpsBlastRunner(outputFolder, this.eval);
        
        DecimalFormat df = new DecimalFormat("#.00");
		df.setRoundingMode(RoundingMode.CEILING);
		
        TSVWriter sig_tsv = new TSVWriter(this.outputFolder+"SpliceVariantResults.tsv", this.build, this.rpsblast);
        
        System.out.println(Utilities.GREEN+"Going through the variants"+ Utilities.RESET);
        int totalVars = vars.size();
        int varsFinished = 0;
        
        while(iter.hasNext()){ //iterate over keys in the vars map
            Map.Entry<String,Variant> entry = iter.next();
            Variant var = entry.getValue();
            var.parseSpliceInfo(rsp, prfr);
            
            System.out.println(var.toString());

          //Run MES and set scores for each variant
            MESRunner mr = new MESRunner(var,this.outputFolder, this.MaxEntPath);  
            if (!mr.IsEmpty()){ // ONLY IF A VALID MES RUN
            	//populate percentDiffList
            	var.checkMesSignificance(); 

            	if (rpsblast.equals("true")) {
            		//Run rpsblast through Cdd Database and find all conserved domains lost
                	rpsRunner.runRPSBlast(var);
            	}
            	
            	// Write out Results
            	sig_tsv.writeVariant(var);
            	
            	
           	} // inside of valid Max Ent Scan only
            
            iter.remove();
            varsFinished = varsFinished + 1;
            
            double progressPercentage = (double) (varsFinished) / (double)(totalVars);
        	final int width = 50; // progress bar width in chars
            
            // update progress
            System.out.print("\r[");
       	    int i = 0;
       	    for (; i < (int)(progressPercentage*width); i++) {
       	      System.out.print(".");
       	    }
       	    for (; i < width; i++) {
       	      System.out.print(" ");
       	    }
       	    System.out.print("] " + df.format(progressPercentage * 100) + "%");
        }
        System.out.print("\n");

        // clean up
        sig_tsv.close();
        Files.deleteIfExists(new File(this.outputFolder+"threePrime.txt").toPath());
        Files.deleteIfExists(new File(this.outputFolder+"fivePrime.txt").toPath());

    }

    private void runAnnotations(String newFile) {
        AnnovarRunner AR = new AnnovarRunner(this.annovar, this.outputFolder, this.build);

        String oneKGenomes = AR.onekGenomes(newFile,this.human);
        GeneralAnnotationParser parser = new GeneralAnnotationParser(this.outputFolder+oneKGenomes);
        this.vars = parser.parse(this.vars);
        
        String exac = AR.Exac(newFile,this.human);
        parser = new GeneralAnnotationParser(this.outputFolder+exac);
        this.vars = parser.parse(this.vars);
        
        String dbscsnv = AR.dbscSNV(newFile,this.human);
        parser = new GeneralAnnotationParser(this.outputFolder+dbscsnv);
        this.vars = parser.parse(this.vars);
        
        if (!this.build.equals("hg38")) {
        	String gerp = AR.Gerp2(newFile, this.human);
        	parser = new GeneralAnnotationParser(this.outputFolder+gerp);
        	this.vars = parser.parse(this.vars);
        }
        
    }

    private String convertAnnovar(){
        AnnovarRunner AR = new AnnovarRunner(this.annovar,this.outputFolder, this.build);
        String avinput = AR.convert2Annovar(this.input);
        String varfunc = AR.Gene(avinput,this.human);
        return varfunc;
    }

    private String writeNewAvinput() {
        try {
            String filename = "SpliceOnly.avinput";
            File newavinput = new File(this.outputFolder + filename);
            FileWriter writer = new FileWriter(newavinput);
            for (Map.Entry<String, Variant> entry : this.vars.entrySet()) {
                String key = entry.getKey();
                String[] chrpos = key.split(":");
                Variant v = entry.getValue();

                writer.write(chrpos[0] + "\t" + chrpos[1] + "\t" + chrpos[1] + "\t" + v.getRef() + "\t" + v.getAlt() + "\n");
            }
            writer.close();

            return filename;

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    private String geneBasedAnnotation(String varfunc){
        annovarWriter nonSplice = new annovarWriter(this.outputFolder+"NonSplice.txt");
        GeneParser gp = new GeneParser(this.outputFolder+varfunc);
        this.vars = gp.parse(nonSplice);
        return writeNewAvinput();
    }
}