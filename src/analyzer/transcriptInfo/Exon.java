package analyzer.transcriptInfo;

/**
 * Created by mwads on 1/26/16.
 */
public class Exon {

    private StringBuilder Seq;
    private String Start;
    private String end;

    public Exon(String Start, String End, StringBuilder Seq, String Strand){
        this.Seq = new StringBuilder(Seq);
        this.Start = Start;
        this.end = End;
    }

    public String getSeq() {
        return Seq.toString();
    }
    public void setSeq(StringBuilder seq) {
        Seq = seq;
    }

    public String getStart() {
        return Start;
    }
    public void setStart(String start) {
        Start = start;
    }

    public String getEnd() {
        return end;
    }
    public void setEnd(String end) {
        this.end = end;
    }


    public Integer getLength(){
        return Seq.length();
    }

}
