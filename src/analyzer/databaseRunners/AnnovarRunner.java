package analyzer.databaseRunners;

import analyzer.Utilities.Utilities;
import java.io.IOException;
/**
 * Created by mwads on 1/15/16.
 */
public class AnnovarRunner {

    private String AnnovarPath;
    private String OutputFolder;
    private String build;

    public AnnovarRunner(String Path2Annovar, String OutputFolder, String build){
        this.AnnovarPath = Path2Annovar;
        this.OutputFolder = OutputFolder;
        this.build = build;
    }

    public String convert2Annovar(String vcf) {
        System.out.println(Utilities.GREEN+"Running VCF to Annovar Conversion"+ Utilities.RESET);
        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/convert2annovar.pl","--splicing_threshold", "35", "-format","vcf4",vcf,"--outfile",vcf+".avinput"};
            ProcessBuilder pb = new ProcessBuilder(call);

            Process p = pb.start();
            p.waitFor();
            System.out.println(Utilities.getProcessOutput(p));
            System.out.println(Utilities.getProcessError(p));
            p.destroy();

            String[] temp = vcf.split("/");

            return temp[temp.length-1]+".avinput";
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;
    }

    public String Gene(String avinput, String human){
        System.out.println(Utilities.GREEN+"Running Gene Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","--splicing_threshold","35","--buildver", build,"-hgvs","-out",this.OutputFolder+avinput,avinput,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            //System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            return avinput+".variant_function";

        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;

    }

    public String onekGenomes(String input, String human){
        System.out.println(Utilities.GREEN+"Running 1000 Genomes MAF Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","-filter","-dbtype","1000g2015aug_all","-buildver", build,"-out",this.OutputFolder+input,this.OutputFolder+input,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            //System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            StringBuilder sb = new StringBuilder(".");
            sb.append(build + "_ALL.sites.2015_08_");
            return input + sb.toString();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;
    }

    public String Gerp2(String input, String human){
        System.out.println(Utilities.GREEN+"Running Gerp++ Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","-filter","-dbtype","gerp++gt2","-buildver", "hg19","-out",this.OutputFolder+input,this.OutputFolder+input,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            //System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            StringBuilder sb = new StringBuilder(".");
            sb.append(build + "_gerp++gt2_");
            return input + sb.toString();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;

    }

    public String Exac(String input, String human){
        System.out.println(Utilities.GREEN+"Running Exac Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","-filter","-dbtype","exac03","-buildver", build,"-out",this.OutputFolder+input,this.OutputFolder+input,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            //System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            StringBuilder sb = new StringBuilder(".");
            sb.append(build + "_exac03_");
            return input + sb.toString();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;

    }
    
    public String dbscSNV(String input, String human){
        System.out.println(Utilities.GREEN+"Running dbscSNV Annotation"+ Utilities.RESET);

        try {
            String[] call = new String[]{"perl",this.AnnovarPath+"/annotate_variation.pl","-filter","-dbtype","dbscsnv11","-buildver", build,"-out",this.OutputFolder+input,this.OutputFolder+input,human};
            ProcessBuilder pb = new ProcessBuilder(call);
            Process p = pb.start();
            p.waitFor();
            //System.out.println(Utilities.getProcessError(p));
            Utilities.getProcessOutput(p);
            p.destroyForcibly();
            p.waitFor();
            StringBuilder sb = new StringBuilder(".");
            sb.append(build + "_dbscsnv11_");
            return input + sb.toString();
        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        return null;

    }
}
