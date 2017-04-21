package analyzer.transcriptInfo;

import analyzer.RefSeq.PullRegionsFromRef;
import analyzer.RefSeq.RefSeqParser;
import analyzer.Utilities.Utilities;

import java.util.*;

/**
 * This class holds the information for the CDS region in question, extracts the CDS sequences, and pulls the sequences
 * used in MaxEntScan.
 *
 * Created by mwads on 1/26/16.
 */
public class CDS {

    private ArrayList<Exon> Exons;
    private ArrayList<Intron> Introns;
    private String Protein;
    //private String geneName;
    private String transName; //'NM_0011'
    private String cdsStart; // column 6 refseq
    private String cdsEnd; // column 7 refseq
    private String chr; // column 2 refseq
    private String strand; // column 3 refseq
    private String[] exonStarts; // column 9 refseq
    private String[] exonEnds; // column 10 refseq
    private String[] exonFrames; // column 15 refseq
    private StringBuilder seq;
    private String cDot; // .vcf.avinput.variant_functions
    private ArrayList<String> CDotList; // ["4451","+","37","C","A"] .vcf.avinput.variant_functions
    private String exon; // .vcf.avinput.variant_functions
    private Integer prime; //this will be either 3 or 5 depending on which side of the intron it is
    private Integer exonSpliceMissed; // the number of the exon that is closest to the offending variant
    private String modifiedProtein;

    public CDS(String transName, String gene) {
        //this.geneName = gene;
        this.transName = transName; //'NM_0011'
        this.Exons = new ArrayList<>();
        this.Introns = new ArrayList<>();
        this.seq = new StringBuilder();


    }

    //Getters and Setters
    public String getStrand() {
        return strand;
    }
    public void setStrand(String strand)
    {
        this.strand = strand;
    }

    public String getcDot() {
        return cDot;
    }
    public void setcDot(String cDot) {
        this.cDot = cDot;
    }

    public ArrayList<String> getCDotList() {
        return CDotList;
    }
    public void setCDotList(ArrayList<String> CDotList) {
        this.CDotList = CDotList;
    }

    public String getExon() {
        return exon;
    }
    public void setExon(String exon) {
        this.exon = exon;
    }

    public String getTransName(){
        return this.transName;
    }
    
    public String getOriginalProtein() {
    	return Protein;
    }
    
    public String getModifiedProtein() {
    	return modifiedProtein;
    }

    // made by CS 4/7/17
    @Override
    public String toString()
    {
    	
/*    private String[] exonStarts; // column 9 refseq
    private String[] exonEnds; // column 10 refseq
    private String[] exonFrames; // column 15 refseq
    private ArrayList<String> CDotList; // ["4451","+","37","C","A"] .vcf.avinput.variant_functions
*/

    	return "\tCDS:\n" +
                "Protein='" + Protein + '\'' +
                ", transName='" + transName + '\'' +
                ", geneName='" + cdsStart + '\'' +
                ", cdsStart='" + cdsStart + '\'' +
                ", cdsEnd='" + cdsEnd + '\'' +
                ", chr=" + chr +
                ", strand=" + strand +
                ", seq='" + seq.toString() + '\'' +
                ", cDot='" + cDot + '\'' +
                ", exon=" + exon +
                ", prime=" + prime +
                ", exonSpliceMissed='" + exonSpliceMissed + '\'' +
                ", modifiedProtein='" + modifiedProtein + '\'' +
                '}';
    	
    }
    
    // made by CS 4/20/17
    public void exonStartStop_printString() {
    	System.out.println("\nExon Starts and Stops:");
    	for (int i = 0; i < exonStarts.length; ++i) 
    	{
    		System.out.println(exonStarts[i] + "-" + exonEnds[i]);
    	}
    }
    
    public void extractCDS(RefSeqParser rsp, PullRegionsFromRef prfr) {
    	//System.out.println("---- CDS ---- extractCDS ----");
    	
    	// info:
    	//0) bin 1) name 2) chrom 3) strand 4) txStart 5) txEnd 6) cdsStart 7) cdsEnd
    	//8) exonCount 9) exonStarts 10) exonEnds 11) score 12) name2 13) cdsStartStat 14) cdsEndStat 15) exonFrames
    	
        String info = rsp.getRefSeqData(transName); //'NM_0011'
        	
        String[] data = info.split("\\t");
        this.chr = data[2];
        this.strand = data[3];
        Integer exonNum = Integer.valueOf(data[8]);
        this.cdsStart = data[6];
        this.cdsEnd = data[7];
        this.exonStarts = data[9].split(",");
        this.exonEnds = data[10].split(",");
        this.exonFrames = data[15].split(",");
        this.seq = new StringBuilder();

        boolean firstExon = true;

        if(this.strand.equals("+"))
        {
            extractCDSRegionPosStrand(prfr,exonNum,firstExon); // CDS:extractCDSRegionPosStrand
        }
        else
        {
            extractCDSRegionNegStrand(prfr,exonNum,firstExon);
        }
        
        
    }

    private void extractCDSRegionNegStrand(PullRegionsFromRef prfr, Integer exonNum, boolean firstExon){
    	//System.out.println("---- CDS ---- extractCDSRegionNegStrand ----");
    	//System.out.println("exonNum: " + exonNum);

        for(int i=exonNum - 1;i >= 0;i--){
            if(Integer.valueOf(this.exonFrames[i]) != -1) { // in the UTR
                if(!firstExon && Integer.valueOf(this.cdsStart) < Integer.valueOf(this.exonStarts[i])){
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.exonStarts[i]), Integer.valueOf(this.exonEnds[i]), this.strand);
                    Exon e = new Exon(this.exonEnds[i],this.exonStarts[i],exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
                }
                else if(firstExon){
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.exonStarts[i]), Integer.valueOf(cdsEnd), this.strand);
                    Exon e = new Exon(this.cdsEnd,this.exonStarts[i],exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
                    firstExon = false;
                }
                else{
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.cdsStart), Integer.valueOf(this.exonEnds[i]), this.strand);
                    Exon e = new Exon(this.exonEnds[i],this.cdsStart,exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
                }

                if(i > 0) {
                    StringBuilder intron = prfr.getRegion(this.chr, Integer.valueOf(this.exonEnds[i - 1]), Integer.valueOf(this.exonStarts[i]), this.strand);
                    Intron I = new Intron(intron, this.exonStarts[i], this.exonEnds[i - 1],this.strand);

                    this.Introns.add(I);
                }
                //System.out.println("\nexonNumber: " + i);
            	//System.out.println("exonStart: " + this.exonStarts[i]);
            	//System.out.println("exonEnd: " + this.exonEnds[i]);

            	//System.out.println("exonLength: " + Exons.get(i).getLength());
            }
        }
    	System.out.println("exonSize: " + Exons.size());

        this.Protein = Utilities.translateProt(this.seq);
    }

    private void extractCDSRegionPosStrand(PullRegionsFromRef prfr, Integer exonNum, boolean firstExon){
    	//System.out.println("---- CDS ---- extractCDSRegionPosStrand ----");
    	//System.out.println("exonNum: " + exonNum);

        for(int i=0;i < exonNum;i++){
            if(Integer.valueOf(this.exonFrames[i]) != -1) { // not in the UTR
                if(!firstExon && Integer.valueOf(this.cdsEnd) > Integer.valueOf(this.exonEnds[i])){
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.exonStarts[i]), Integer.valueOf(this.exonEnds[i]), this.strand);
                    Exon e = new Exon(this.exonStarts[i],this.exonEnds[i],exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
                }
                else if(firstExon){
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(cdsStart), Integer.valueOf(this.exonEnds[i]), this.strand);
                    Exon e = new Exon(this.cdsStart,this.exonEnds[i],exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());

                    firstExon = false;
                }
                else{
                    StringBuilder exon = prfr.getRegion(this.chr, Integer.valueOf(this.exonStarts[i]), Integer.valueOf(this.cdsEnd), this.strand);
                    Exon e = new Exon(this.exonStarts[i],this.cdsEnd,exon,this.strand);
                    this.Exons.add(e);
                    this.seq.append(exon.toString());
                }

                if(i < exonNum-1) {
                    StringBuilder intron = prfr.getRegion(this.chr, Integer.valueOf(this.exonEnds[i]), Integer.valueOf(this.exonStarts[i + 1]), this.strand);
                    Intron I = new Intron(intron, this.exonEnds[i], this.exonStarts[i + 1],this.strand);

                    this.Introns.add(I);
                }
                
            	//System.out.println("\nexonNumber: " + i);
            	//System.out.println("exonStart: " + this.exonStarts[i]);
            	//System.out.println("exonEnd: " + this.exonEnds[i]);

            	//System.out.println("exonLength: " + Exons.get(i).getLength());
            	
            }
        }
    	System.out.println("exonSize: " + Exons.size());

        this.Protein = Utilities.translateProt(this.seq);
    }

    public String getMES3Prime(Integer pos){ // creates 23 base sequence (20 intron + 3 exon)
    	
            StringBuilder sb = new StringBuilder();
            Integer e = -1;
            e = getPosExon(pos);

            this.exonSpliceMissed = e;
            this.prime = 3;

            StringBuilder intron = new StringBuilder(this.Introns.get(e-1).getSeq());
            sb.append(intron.substring(intron.length()-20,intron.length()));
            StringBuilder seq = new StringBuilder(this.Exons.get(e).getSeq());
            if (seq.length() < 3) {
            	sb = new StringBuilder("SEQ TOO SHORT");
            } else {
                sb.append(seq.substring(0,3));
            }
            return sb.toString();
    }
    
    public String getMES5Prime(Integer pos){ // creates 9 base sequence (3 exon + 6 intron)
    	
            StringBuilder sb = new StringBuilder();
            Integer e = -1;
            e = getPosExon(pos);

            this.exonSpliceMissed = e;
            this.prime = 5;
            StringBuilder seq = new StringBuilder(this.Exons.get(e).getSeq());
            sb.append(seq.substring(seq.length()-3,seq.length()));
            StringBuilder intron = new StringBuilder(this.Introns.get(e).getSeq());
            sb.append(intron.substring(0,6));
            if (intron.length() < 6) {
            	sb = new StringBuilder("SEQ TOO SHORT");
            } else {
                sb.append(intron.substring(0,6));
            }
            return sb.toString();

    }

    private Integer getPosExon(Integer Pos){ // get the position of the exon in the sequence
    	//System.out.println("---- CDS ---- getPosExon ----");

        Integer total = 0;
        Integer j;
        Integer curr_exon = -1;
        System.out.println("position: " + Pos);

        for(j = 0; j < this.Exons.size(); ++j){
            total += this.Exons.get((j)).getLength();
            System.out.println("total: " + total);
            System.out.println("j: " + j);

            if(Pos <= total) {
                return j;
            } else {
            	curr_exon = j;
            }
            	
        }
        System.out.println("exon number=" + j);
        System.out.println("curr_exon  =" + curr_exon);

        System.out.println("  exon size=" + this.Exons.size());

        return curr_exon;
    }


    public String makeModifiedProtein(){
    	//System.out.println("---- CDS ---- makeModifiedProtein ----");

        StringBuilder sb = new StringBuilder();
        for(int i = 0; i < Exons.size(); i++){
            if(i == exonSpliceMissed && prime == 3){
                sb.append(Introns.get(exonSpliceMissed-1).getSeq());
            }
            sb.append(Exons.get(i).getSeq());
            if(i == exonSpliceMissed && prime == 5){
                sb.append(Introns.get(exonSpliceMissed).getSeq());
            }
        }
        System.out.println("seq: " + sb.toString());
        modifiedProtein = Utilities.translateProt(sb);

        return modifiedProtein;
    }


}
