package analyzer;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentChoice;
//import net.sourceforge.argparse4j.inf.ArgumentChoice;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

/* 
 * This class is used to create an argument parser for command line arguments to passed in,
 *  and to run the algorithm. 
 * It contains the main method. 
*/
public class SpliceEngine {

    public static void main( String[] args)
    {
        SpliceEngine se = new SpliceEngine();
        ArgumentParser parser = se.setup_parser(); //Set up the argument parser 
        try {
            Namespace res = parser.parseArgs(args);
            SpliceRunner sr = new SpliceRunner(res); 
            sr.run(); // Run the algorithm with the arguments passed in. 
        } catch (ArgumentParserException e) {
            parser.handleError(e);
        } catch (Exception e){
            e.printStackTrace();
        }
    }


    //-------------------------------------------------------------------------------------
    // This method used to create a command line argument parser, which will set up required 
    //  arguments and flags to be entered through running the the algorithm on the command line
    // Flags:
    // -i = the input file in VCF format 
    // -o = the name of the directory to save the output to
    // -a = path to annovar 
    // -d = path to the HumandDB (databases like grep++, exac, and 1kgenomes) downloaded using annovar. 
    // -h = path to the Reference Genome from USCS (ex = hg19)
    // -r = path to file that contains refSeq data from UCSC table browser 
    // -m = the full path to the MaxEntScan directory
    // -s = path to Samtools // assumed to be in path
    //-------------------------------------------------------------------------------------
    private ArgumentParser setup_parser(){
        ArgumentParser parser = ArgumentParsers.newArgumentParser("SpliceSiteIdentifier.py")
                .description("This program analyzes a vcf file to find and score Splice Site Variants");

        //REQUIRED
        
        parser.addArgument("-i","--input")
  	  		  .dest("Input")
  	  		  .help("This is the input file (VCF format only).")
  	  		  .required(true)
  	  		  .type(String.class);

        parser.addArgument("-o","--output")
   	  		  .dest("Output")
   	  		  .help("This is the the name of the directory to which the output files are saved.")
   	  		  .required(true)
   	  		  .type(String.class); 
  
        parser.addArgument("-a","--annovar")
        	  .dest("Annovar")
        	  .help("This is the path to annovar.")
        	  .required(true)
        	  .type(String.class);

        parser.addArgument("-d","--humandb")
         	  .dest("human")
         	  .help("This is the path to the humandb that Annovar uses.")
              .required(true)
              .type(String.class);
        
        parser.addArgument("-m","--MaxEntScan")
        	  .dest("MaxEntPath")
        	  .help("This is the full path to the MaxEntScan directory.")
        	  .required(true)
        	  .type(String.class);
        
        parser.addArgument("-g","--genome")
          	  .dest("Ref")
          	  .help("This is the path to the directory that contains the UCSC reference genome by chromosome downloaded (hg19/hg38).")
          	  .required(true)
          	  .type(String.class);
      

        // NOT REQUIRED (DEFAULTS SET)
        parser.addArgument("-b","--buildver")
        	  .dest("build")
        	  .required(false)
  	  		  .setDefault("hg19")
  	  		  .type(String.class)
        	  .choices(new ArgumentChoice() {
        		  @Override
        		  public boolean contains(Object o) {
        			  String val = String.valueOf(o);
        			  if (val.equals("hg19")) {
        				  return true;
        			  } else if (val.equals("hg38")) {
        				  return true;
        			  }
        			  return false;
        		  }

        		  @Override
        		  public String textualFormat() {
        			  return "Valid options are 'hg19'(default) or 'hg38'";
        		  }
        	  });
        	  
        parser.addArgument("-r","--RefSeqFile")
        	  .dest("analyzer/RefSeq")
        	  .help("This is the path to the file that contains the UCSC table viewer RefSeq data.")
        	  .required(false)
        	  .setDefault("RefSeq_hg19.txt")
        	  .type(String.class);  
        
        parser.addArgument("-e","--evalue")
  	  		  .dest("eval")
  	  		  .help("This is the evalue cut off for the rpsblast of the Cdd database.")
  	  		  .required(false)
  	  		  .setDefault(".005")
  	  		  .type(String.class);  
        
        parser.addArgument("-s","--Samtools")
         	  .dest("Samtools")
         	  .help("This is the path to the Samtools executable.")
         	  .required(false)
         	  .setDefault("samtools") // in path assumed
         	  .type(String.class);
        
        parser.addArgument("-u","--debug")
   	  		  .dest("debug")
   	  		  .help("Sets debug mode")
   	  		  .required(false)
   	  		  .setDefault("false") // in path assumed
   	  		  .type(String.class);


       
        return parser;
    }
}


/*parser.addArgument("-R","--spliceSiteRange")
.dest("range")
.help("Give the number of bases that you want defined as a splice site. It will look that many bases to either side of the splice site.")
.type(Integer.class);
*/
/*parser.addArgument("-S","--scoringAlgorithm").dest("score").choices(new ArgumentChoice() {
@Override
public boolean contains(Object o) {
String val = String.valueOf(o);
if (val.equals("MaxEntScan") || val.equals("MES")) {
  SpliceEngine.algorithm = "MES";
  return true;
} else if (val.equals("Ensemble") || val.equals("EN")) {
  SpliceEngine.algorithm = "EN";
  return true;
}
return false;
}

@Override
public String textualFormat() {
return "The choices are \"MaxEntScan\" aka \"MES\", and \"Ensemble\" aka \"EN\"";
}
}).required(true).help("Choose the desired algorithm.");
*/


// ************ Example Bash script to run the algorithm *********************
//java -jar \
//         /user/splice-variant-analyzer/target/SVA-1.0-SNAPSHOT-jar-with-dependencies.jar \
//         -S MES \
//         -A ~/user/software/annovar \
//         -H ~/user/compute/humandb/ \
//         -F ~/user/hg19Ref/ \ 
//         -R 50 \
//         -i /user/NA12878.chr21.vcf \
//         -o outfolder/ \
//         -s ~/user/software/samtools-1.0/samtools \
//         -D ~/user/refSeqData \
//         -a /user/maxent/
//***************************************************************************

//   CORINNE'S!
//************ Example Bash script to run the algorithm *********************
//java -jar \
//      /user/splice-variant-analyzer/target/SVA-1.0-SNAPSHOT-jar-with-dependencies.jar \
//      -S MES \
//      -A /Users/coripenrod/software/annovar \
//      -H /fslhome/RidgeLab/compute/humandb/ \
//      -F /Users/coripenrod/hg19/ \ 
//      -R 50 \
//      -i /Users/coripenrod/software/new_test.vcf \
//      -o /Users/coripenrod/COLLEGE/RidgeLab/SVA/ \
//      -s /Users/coripenrod/software/samtools-1.4/samtools \
//      -D /fslhome/penrodce/refSeqData \comes from UCSC genome table viewer (assembly = hg19, track = refSeq genes, table = refGene, cdsstart and end not noce(filter))
//      -a /Users/coripenrod/maxent/
//***************************************************************************
