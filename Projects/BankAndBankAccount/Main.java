package Exercise3;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.TreeSet;

public class Main {

    public static void main(String[] args) {


        Bank bcr = new Bank("Banca Comerciala Romana");

        for (int i = 0; i < 20; ++i){
            String myString = "Razvan" + i;

            bcr.addAccount(myString, (int) (Math.random() * 100));
        }
        System.out.println("<--------------------SORTED LIST BY AMOUNT-------------------------->");
        bcr.printAccounts();

        System.out.println("<--------------------SORTED LIST BETWEEN AMOUNTS-------------------------->");

        bcr.printAccounts(20,40);

        System.out.println("<--------------------SORTED LIST ALPHABETICALLY-------------------------->");

        sortAlphabetically(bcr.getMembers());

    }



    public static void sortAlphabetically(TreeSet<BankAccount> myInitialSet) {

        TreeSet<BankAccount> sortedSet = new TreeSet<>(new Comparator<BankAccount>() {
            @Override
            public int compare(BankAccount o1, BankAccount o2) {
                return o1.getOwner().compareTo(o2.getOwner());
            }
        });

        Iterator iterator = myInitialSet.iterator();
        while (iterator.hasNext()){
            BankAccount b = (BankAccount) iterator.next();
           sortedSet.add(b);
        }



        iterator = sortedSet.iterator();
        while (iterator.hasNext()){
            BankAccount b1 = (BankAccount) iterator.next();
            System.out.println(b1.toString());
        }


    }
}
