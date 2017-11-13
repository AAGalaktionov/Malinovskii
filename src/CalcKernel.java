import com.wolfram.jlink.KernelLink;
import com.wolfram.jlink.MathLinkFactory;

/**
 * Created by Alex on 05.11.2017.
 */

public class CalcKernel {
    static KernelLink ml;
    static String path = "-linkmode launch -linkname 'C:/Program Files/Wolfram Research/Mathematica/11.2/MathKernel'";

    public CalcKernel() {


        try {
            ml = MathLinkFactory.createKernelLink(path);// подключаем ядро
            ml.discardAnswer();// дожидаемся загрузки ядра
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public double calculate(String s) {

        try {
            String toCalculate = "N[%s, 50]";
            ml.evaluate(String.format(toCalculate, s));
            ml.waitForAnswer();
            return ml.getDouble();// считываем результат;
        } catch (Exception e) {
        e.printStackTrace();
        throw new RuntimeException(e);
        }
    }
}
