package analyzer.variantInfo;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import analyzer.RefSeq.PullRegionsFromRef;
import analyzer.RefSeq.RefSeqParser;
//import analyzer.Utilities.Utilities;
//import analyzer.fileWriters.VCFWriter;
import analyzer.transcriptInfo.CDS;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Created by mwads on 1/26/16.
 */
public class Variant {

    private String Chr; // chr21
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
    private ArrayList<Double> percentDiffList;
    public ArrayList<Integer> WithinGenePos;    // from cDot list
    
    public ArrayList<String> ConservedDomains;
    public ArrayList<String> PDBList;


    public Variant(String chr, Integer pos, String spliceInfo, String ref, String alt, String homhet){
        this.Chr = chr;
        this.Pos = pos;
        this.Ref = ref;
        this.Alt = alt;
        this.homhet = homhet;
        this.spliceInfo = spliceInfo;
        this.Annotations = new ArrayList<>();
        this.transcripts = new ArrayList<>();
        this.OriginalMesScores = new ArrayList<>();
        this.VariantMesScores = new ArrayList<>();
        this.percentDiffList = new ArrayList<>();
        this.WithinGenePos = new ArrayList<>();


        this.ConservedDomains = new ArrayList<>();
        this.PDBList = new ArrayList<>();


    }

   // Getters and Setters

    public String getChr() {
        return Chr;
    }
    public void setChr(String chr) {
        Chr = chr;
    }

    public Integer getPos() {
        return Pos;
    }
    public void setPos(Integer pos) {
        Pos = pos;
    }

    public String getRef() {
        return Ref;
    }
    public void setRef(String ref) {
        Ref = ref;
    }

    public String getAlt() {
        return Alt;
    }
    public void setAlt(String alt) {
        Alt = alt;
    }

    public String getGeneName() {
        return GeneName;
    }
    
    public String getHomhet() {
        return homhet;
    }
    public void setHomhet(String homhet) {
        this.homhet = homhet;
    }

    public String getSpliceInfo() {
        return spliceInfo;
    }
    public void setSpliceInfo(String spliceInfo) {
        this.spliceInfo = spliceInfo;
    }

    public ArrayList<Double> getOriginalMesScores() {
        return OriginalMesScores;
    }
    public void setOriginalMesScores(ArrayList<Double> originalMesScores) {
        OriginalMesScores = originalMesScores;
    }
    
    public ArrayList<Double> getVariantMesScores() {
        return VariantMesScores;
    }
    public void setVariantMesScores(ArrayList<Double> variantMesScores) {
        VariantMesScores = variantMesScores;
    }
    
    public ArrayList<CDS> getCDSList() {
    	return transcripts;
    }
    
	public ArrayList<Double> getPercentDiffList() {
		return percentDiffList;
	}
    
    // Creating variants to write

    public VariantContextBuilder createVariantContext(){
    	//System.out.println("---- Variant ---- createVariantContext ----");

        ArrayList<Allele> alleles = new ArrayList<>();
        alleles.add(Allele.create(this.Ref,true));
        alleles.add(Allele.create(this.Alt,false));
        VariantContextBuilder builder = new VariantContextBuilder("MES Filtering", this.Chr, Long.valueOf(this.Pos), Long.valueOf(this.Pos), alleles);

        builder.attribute("Gene",this.GeneName);

        builder.noGenotypes();
        return builder;
    }

    // Gathering info for variant

    public void addAnnotation(String annotation){
        this.Annotations.add(annotation);
    }

    @Override
    public String toString() {
        /*return "Variant{" +
                "Chr='" + Chr + '\'' +
                ",\nPos=" + Pos +
                ",\nRef='" + Ref + '\'' +
                ",\nAlt='" + Alt + '\'' +
                ",\nhomhet='" + homhet + '\'' +
                ",\nspliceInfo='" + spliceInfo + '\'' +
                ",\nAnnotations=" + String.valueOf(Annotations.size()) +
                ",\ntranscripts=" + String.valueOf(transcripts.size()) +
                ",\nGeneName='" + GeneName + '\'' +
                ",\nVariantMESScores='" + VariantMES
                "}\n\n";
                */
    	StringBuilder sb = new StringBuilder();
    	for (CDS cd : transcripts) {
    		sb.append(cd.toString());
    	}
    	
    	StringBuilder sb2 = new StringBuilder();

    	for (Double d : VariantMesScores) {
    		sb2.append(d.toString());
    	}
    	
    	return "Variant{\n" +
        "\tChr='" + Chr + '\'' +
        ",\tPos=" + Pos + ",\n" +
        "\tspliceInfo='" + spliceInfo + "\',\n" +
        "\ttranscripts=" + String.valueOf(transcripts.size()) +
        "\tMesScores=" + String.valueOf(VariantMesScores.size()) +
        //",\n" + sb.toString() +
        ",\n" + sb2.toString() +
        "}\n\n";
    }

    public void parseSpliceInfo(RefSeqParser rsp, PullRegionsFromRef prfr){ // info from the varfunct (.vcf.avinput.variant_function)
    	
        String[] geneList = this.spliceInfo.split("\\),"); // DIP2A(NM_001146116:exon37:c.4451+37C>A,NM_015151:exon37:c.4463+37C>A)
        for(String gene : geneList) { // ['GENE', 'NM_0011:exon47:c.4451+37C>A,NM_015151:exon37:c.4463+37C>A']
            String[] gene2trans = gene.split("\\("); // ['GENE', 'NM_0011:exon47:c.4451+37C>A,NM_015151:exon37:c.4463+37C>A']
            this.GeneName = gene2trans[0];
            if (gene2trans.length == 1)
                return;

            String[] transInfo = gene2trans[1].replace(")","").split(","); // ['NM_0011:exon47:c.4451+37C>A', 'NM_015151:exon37:c.4463+37C>A']
            for (String trans : transInfo) { // NM_0011:exon47:c.4451+37C>A,NM_015151:exon37:c.4463+37C>A
                String[] info = trans.split(":"); // ['NM_0011', 'exon47', 'c.4451+37C>A']
                if (rsp.getRefSeqData(info[0]) != null) {; // if it is in the refSeq Data!
                	CDS cds = new CDS(info[0], gene2trans[0]);
                	cds.setcDot(info[2]); // c.4451+37C>A
                	cds.setCDotList(parseCDot(info[2])); // ["4451","+","37","C","A"]
                	
			if (cds.getCDotList().size() == 4) {
                		if (cds.getCDotList().get(1).equals("+")) {
                    		this.WithinGenePos.add(Integer.parseInt(cds.getCDotList().get(0)) + Integer.parseInt(cds.getCDotList().get(2)));
                		} else {
                			this.WithinGenePos.add(Integer.parseInt(cds.getCDotList().get(0)) - Integer.parseInt(cds.getCDotList().get(2)));
                		}		
                				
                		cds.setExon(info[1]);
                		cds.extractCDS(rsp, prfr);
                		
                		// somehow check for duplicate transcripts?
                		
                		/*boolean dupl = false;
                		for (CDS t : transcripts){
                			System.out.println(t.cDot + ": " + t.seq);
                			System.out.println(cds.cDot + ": " + cds.seq);

                			if (t.seq.equals(cds.seq)) {
                				System.out.println("\nDUPLICATE\n");
                				t.addDupl(cds);
                				dupl = true;
                			}
                		}*/
                		//if (!dupl) 
                		this.transcripts.add(cds);             	
		    }
                }
            }
        }
    }

    private ArrayList<String> parseCDot(String cDot){
        Pattern p = Pattern.compile("c\\.(\\d+)([-+])(\\d+)([ATGC])>([ACTG])");
        Matcher m = p.matcher(cDot);
        ArrayList<String> info = new ArrayList<>();
        if(m.find()){
            info.add(m.group(1));
            info.add(m.group(2));
            info.add(m.group(3));
            info.add(m.group(4));
            info.add(m.group(5));
        }
        return info;
    }

    public String makeMESSequenceFile(String outFolder) throws Exception{
    	
        VariantContextBuilder vcb = createVariantContext();

        String threePrime = null;
        String fivePrime = null;
        boolean three = false;
        boolean five = false;

        FileWriter threePrimeFile = null;
        FileWriter fivePrimeFile = null;

        try {
            threePrimeFile = new FileWriter(outFolder + "threePrime.txt");
            fivePrimeFile = new FileWriter(outFolder + "fivePrime.txt");
        } catch (IOException e) {
            e.printStackTrace();
        }
        

        LinkedList<String> filteredNames = new LinkedList<>();
        
        for (CDS cds : transcripts){

            List<String> CDotList = cds.getCDotList(); // ["4451","+","11","C","A"]
            if(CDotList.get(1).equals("-")){
                if ( Integer.valueOf(CDotList.get(2)) <= 20) {
                    threePrime = new String(outFolder+"threePrime.txt");
                    three = write3Prime(threePrimeFile,cds);
                }
                else{
                    filteredNames.add(cds.getTransName());
                }
            }
            else
            {
                if(Integer.valueOf(CDotList.get(2)) <= 6) {
                    fivePrime = new String(outFolder+"fivePrime.txt");  
                    five = write5Prime(fivePrimeFile,cds);
                }
                else{
                    filteredNames.add(cds.getTransName());
                }
            }
        }


        if(filteredNames.size()!=0) {
            StringBuilder sb = new StringBuilder();
            for (String t: filteredNames){
                sb.append(t+",");
            }
            sb.replace(sb.length()-1,sb.length(),"");
            vcb.attribute("Transcripts", sb.toString());
        }

        if(threePrime!= null && fivePrime != null)
            throw new Exception("The variant shouldn't end up on both ends of the intron.");
        
        if(!three) {
            if(five) {
                try {
                    fivePrimeFile.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
                return fivePrime;
            }
            else{
                return null;
            }
        }
        else if (!five){
            try {
                threePrimeFile.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
            return threePrime;

        } else {
            return null;
        }
        
        
    }

    private boolean write3Prime(FileWriter threePrime, CDS cds){
    	
        try {
            List<String> CDotList = cds.getCDotList();
            String originalSeq = cds.getMES3Prime(Integer.valueOf(CDotList.get(0)));
            StringBuilder sb = new StringBuilder(originalSeq);
            
            if (originalSeq.length() != 13) {
            	sb.setCharAt(20-Integer.valueOf(CDotList.get(2)), CDotList.get(4).charAt(0));
                threePrime.write(">original\n" + originalSeq + "\n>newSeq\n" + sb.toString() + "\n");
                return true;
            } else {
            	return false;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return false;
    }

    private boolean write5Prime(FileWriter fivePrime, CDS cds){
    	
        try {
            List<String> CDotList = cds.getCDotList();
            String originalSeq = cds.getMES5Prime(Integer.valueOf(CDotList.get(0)));
            StringBuilder sb = new StringBuilder(originalSeq);

            if (originalSeq.length() != 13) {
            	sb.setCharAt(Integer.valueOf(2 + Integer.valueOf(CDotList.get(2))), CDotList.get(4).charAt(0));
                fivePrime.write(">original\n" + originalSeq + "\n>newSeq\n" + sb.toString() + "\n");
                return true;
            } else {
            	return false;
            }
            
        } catch (IOException e) {
            e.printStackTrace();
        }
		return false;
    }

    public void makeModifiedProtein(){
        for(CDS cds : transcripts){
            cds.makeModifiedProtein();
        }
    }

    public int checkMesSignificance(){
    	
        this.percentDiffList = new ArrayList<>();
        int sigCount = 0;
        int likelySigCount = 0;
        int notSigCount = 0;
        for(int i=0; i < OriginalMesScores.size(); i++){
            Double percentDiff = ((Double.valueOf(VariantMesScores.get(i))-Double.valueOf(OriginalMesScores.get(i))) / Math.abs(Double.valueOf(OriginalMesScores.get(i))) * 100);
            percentDiffList.add(percentDiff);
            
           if(percentDiff <= -20)
                sigCount++;
            else if(percentDiff <= -10)
                likelySigCount++;
            else
                notSigCount++;
        }
        
        if(sigCount > (likelySigCount+notSigCount))
            //if(this.Annotations.get(0).equals("NA")&&this.Annotations.get(1).equals("NA")&&this.Annotations.get(2).equals("NA")) //This checks to make sure that the variant doesn't have any annotations from 1KG, Exac, or Gerp++
                return 2;
        if(sigCount == (likelySigCount+notSigCount) || likelySigCount > (sigCount+notSigCount))
            return 1;
        return 0;
    }


}
