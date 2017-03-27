package analyzer;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.inf.ArgumentChoice;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;

/* 
 * This class is used to create an argument parser for command line arguments to passed in,
 *  and to run the algorithm. 
 * It contains the main method. 
*/
public class SpliceEngine {

    private static String algorithm;

    public static void main( String[] args)
    {
        SpliceEngine se = new SpliceEngine();
        ArgumentParser parser = se.setup_parser(); //Set up the argument parser 
        try {
            Namespace res = parser.parseArgs(args);
            SpliceRunner sr = new SpliceRunner(res, SpliceEngine.algorithm); 
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
    // -S = which scoring algorithm to use. (Currently only supports MES) 
    // -A = path to annovar 
    // -H = path to the HumandDB created using annovar. 
    // -F = path to the Reference Genome from USCS (ex = hg19)
    // -D = path to file that contains analyzer.refseq data 
    // -S = path to Samtools 
    // -R = define the splice site by the number of bases away from an exon-intron region 
    // -i = the input file in VCF format 
    // -o = the name of the directory to save the output to
    // -a = the full path to the algorithm directory
    //-------------------------------------------------------------------------------------
    private ArgumentParser setup_parser(){
        ArgumentParser parser = ArgumentParsers.newArgumentParser("SpliceSiteIdentifier.py")
                .description("This program analyzes a vcf file to find and score Splice Site Variants");

        parser.addArgument("-S","--scoringAlgorithm").dest("score").choices(new ArgumentChoice() {
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
                return "The choices are \"MaxEntScan\" aka \"MES\", and \"Ensemble\" aka \"EN\")";
            }
        }).required(true).help("Choose the desired algorithm.");



        parser.addArgument("-A","--Annovar").dest("Annovar").help("This is the path to annovar.").required(true).type(String.class);

        parser.addArgument("-H","--humandb").dest("human").help("This is the path to the humandb that Annovar uses.")
                .required(true).type(String.class);

        parser.addArgument("-F","--RefFile").required(true).dest("Ref").help("This is the path to the directory" +
                " that contains the UCSC reference genome by chromosome downloaded.").type(String.class);

        parser.addArgument("-D","--RefSeqFile").required(true).dest("analyzer/RefSeq").help("This is the path to the file" +
                " that contains the analyzer.RefSeq data.").type(String.class);

        parser.addArgument("-s","--Samtools").required(true).dest("Samtools").help("This is the path to the Samtools" +
                "executable.").type(String.class);

        parser.addArgument("-R","--spliceSiteRange").dest("range").help("Give the number of bases that you want " +
                "defined as a splice site. It will look that many bases to either side of the splice site.").type(Integer.class);

        parser.addArgument("-i","--input").dest("Input").type(String.class).help("Here you provide the input in VCF format only.").required(true);

        parser.addArgument("-o","--output").dest("Output").type(String.class).help("Here you provide the name of the directory to which I will save the output files.").required(true);

        parser.addArgument("-a","--algorithm").dest("AlgorithmPath").type(String.class).help("This is the full path to the algorithm directory.").required(true);

        return parser;
    }
}


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
//      		-H /coripenrod/compute/humandb/ \
//      -F /Users/coripenrod/hg19/ \ 
//      -R 50 \
//      -i /Users/coripenrod/software/new_test.vcf \
//      -o /Users/coripenrod/COLLEGE/RidgeLab/SVA/ \
//      -s /Users/coripenrod/software/samtools-1.4/samtools \
//      		-D /user/refSeqData \
//      -a /Users/coripenrod/maxent/
//***************************************************************************
