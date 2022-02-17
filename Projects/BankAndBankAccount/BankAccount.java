package Exercise3;

public class BankAccount implements Comparable {

    private String owner;
    private double balance;

    public BankAccount(String owner, double balance) {
        this.owner = owner;
        this.balance = balance;
    }

    public String getOwner() {
        return owner;
    }

    public double getBalance() {
        return balance;
    }

    public void setOwner(String owner) {
        this.owner = owner;
    }

    public void setBalance(double balance) {
        this.balance = balance;
    }

    public void withdraw(double amount){
        if (amount >= balance){
            System.out.println("Not enough balance to withdraw this amount");
        } else{
            balance -= amount;
            System.out.println("Amount withdrawn successfully");
        }
    }

    public void deposit(double amount){
        balance += amount;
        System.out.println("Amount added successfully");
    }

    @Override
    public int hashCode() {
        return owner.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null || (!(obj instanceof BankAccount)))
            return false;

        BankAccount bankAccount = (BankAccount) obj;
        return bankAccount.getBalance() == this.balance && bankAccount.getOwner().equals(this.getOwner());
    }

    @Override
    public String toString() {
        return "BankAccount { " +
                "owner = " + owner +
                ", balance = " + balance +
                " }";
    }

    @Override
    public int compareTo(Object o){
        BankAccount b1 = (BankAccount) o;
        if (b1.balance > this. balance) return -1;
        if (b1.balance < this.balance) return  1;
        else return 0;
    }


}
