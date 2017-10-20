import java.lang.reflect.Array;
import java.math.*;
import java.util.ArrayList;

public class Main {


    public static void main(String[] args) {
        System.out.println("Begin");
        double alpha = 1, betta = 1, M, D2, u =10 , c = 0, v = 0, t = 100000;
        double F = 0.0398;
        M = alpha/betta;
        D2 = (2*alpha)/(Math.pow(betta,2));
        double  x1 = 1 ;


        ArrayList <Double> Mt = new ArrayList <Double> ();
        while (c <= 2.05) {
            double x2 = ((c*(t-v))/(u+c*v))+1;
            Mt.add( F*(((Math.sqrt(u+c*v))/(c*Math.sqrt(D2)*Math.sqrt(x2)))*(x2*(1-c*M)-1))+
                    Math.exp(2*((u+c*v)/Math.pow(c,2)*D2)*(1-c*M))*F*(1-(Math.sqrt(u+c*v))/(c*Math.sqrt(D2)*Math.sqrt(x2)))*(x2*(1-c*M))-
                    F*(((Math.sqrt(u+c*v))/(c*Math.sqrt(D2)*Math.sqrt(x1)))*(x1*(1-c*M)-1))+
                    Math.exp(2*((u+c*v)/Math.pow(c,2)*D2)*(1-c*M))*F*(1-(Math.sqrt(u+c*v))/(c*Math.sqrt(D2)*Math.sqrt(x1)))*(x1*(1-c*M)));
            c += 0.05;
        }
        for (int j = 0; j < Mt.size(); j++) {
            System.out.println(Mt.get(j));
        }






    }
}
