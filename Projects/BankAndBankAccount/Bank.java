package Exercise3;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.TreeSet;

public class Bank {
    private String name;
     TreeSet<BankAccount> members;

    public Bank(String name) {
        this.name = name;
        members = new TreeSet<>();
    }

    public void addAccount(String owner, double balance){

        BankAccount newBankAccount = new BankAccount(owner, balance);
        members.add(newBankAccount);
    }

    public void printAccounts(){
        Iterator iterator = members.iterator();
        while (iterator.hasNext()){
            System.out.println(((BankAccount) iterator.next()).toString());
        }
    }

    public void printAccounts(double minBalance, double maxBalance){


        Iterator iterator = members.iterator();
        while (iterator.hasNext()){
            BankAccount b = (BankAccount) iterator.next();
            BankAccount b1 = new BankAccount(b.getOwner(),b.getBalance());
            if (b1.getBalance() >= minBalance && b1.getBalance() <= maxBalance){
                System.out.println(b1.toString());
            }
        }
    }

    public BankAccount getAccount(String owner){

        for (Object obj: members){
            if (((BankAccount)obj).getOwner().equals(owner) )
                return (BankAccount) obj;
        }
        return null;
    }

    public String getName() {
        return name;
    }

    public TreeSet<BankAccount> getMembers() {
        return members;
    }
}
