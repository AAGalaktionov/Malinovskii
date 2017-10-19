import java.util.ArrayList;
import java.util.DoubleSummaryStatistics;

public class Main {


    public static void main(String[] args) {
        System.out.println("Begin");
        double alpha = 1, betta = 1, M, D2, u =10 , c = 0.0, v = 0, t = 100;
        double F = 0.0398;
        M = alpha/betta;
        D2 = (2*alpha)/(Math.pow(betta,2));
        double x1 = 1 ;
        int i = 0;
        ArrayList<Double> Mt = new ArrayList <Double>();


        while (c < 2.1) {
           double x2 = ((c*(t-v))/(u+c*v))+1;
            double tmp = (F*((Math.sqrt(u+c*v))/(c*Math.sqrt(D2)*Math.sqrt(x2)))*(x2*(1-c*M)-1)+
                    Math.exp(2*((u+c*v)/Math.pow(c,2.0)*D2)*(1-c*M))*F*(-(Math.sqrt(u+c*v))/(c*Math.sqrt(D2)*Math.sqrt(x2)))*(x2*(1-c*M)+1))-
                    (F*((Math.sqrt(u+c*v))/(c*Math.sqrt(D2)*Math.sqrt(x1)))*(x1*(1-c*M)-1)+
                    Math.exp(2*((u+c*v)/Math.pow(c,2.0)*D2)*(1-c*M))*F*(-(Math.sqrt(u+c*v))/(c*Math.sqrt(D2)*Math.sqrt(x1)))*(x1*(1-c*M)+1));
            c += 0.05;
            Mt.add(tmp);
            i++;
        }
        for (int j = 0; j < Mt.size() ; j++) {
            System.out.println(Mt.get(j));
        }
                System.out.println("end");
    }
}
