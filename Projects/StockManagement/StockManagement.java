package exercitiul2;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import java.util.*;

public class StockManagement{
	public static void main(String[] args) {
		JFrame management= new Controller();
		management.setVisible(true);
	}
}

class Controller extends JFrame{
	
	public ArrayList<Product> list= new ArrayList<Product>();
	
	private JButton addProduct= new JButton("Add new product");
	private JButton viewProducts= new JButton("View available products");
	private JButton deleteProduct=new JButton ("Delete product");
	private JButton changeProductAvailableQuantity= new JButton("Change product available quantity");
	private JLabel text= new JLabel("Welcome to the stock management!");
	
	Controller(){
		super("Stock Management");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        init();
        setSize(540,400);
        setVisible(true);
	}
	
	public void init() {
        this.setLayout(null);
        text.setBounds(160,15,250, 50);
        add(text);
        addProduct.setBounds(150, 70, 230, 40);
        add(addProduct);
        addProduct.addActionListener(new TreatButtonAddProduct());
        viewProducts.setBounds(150, 130, 230, 40);
        add(viewProducts);
        viewProducts.addActionListener(new TreatButtonViewProducts());
        deleteProduct.setBounds(150, 190, 230, 40);
        add(deleteProduct);
        deleteProduct.addActionListener(new TreatButtonDeleteProduct());
        changeProductAvailableQuantity.setBounds(150, 250, 230, 40);
        add(changeProductAvailableQuantity);
        changeProductAvailableQuantity.addActionListener(new TreatButtonChangeQuantity());
	}
	
	class TreatButtonAddProduct implements ActionListener{
		public void actionPerformed(ActionEvent e) {
			secondFrameAdd s= new secondFrameAdd(list);
			s.setVisible(true);
		}
	}
	
	class TreatButtonViewProducts implements ActionListener{
		public void actionPerformed(ActionEvent e) {
			if(list.isEmpty()==true) { 
				JOptionPane.showMessageDialog(null, 
						"There are no products! Please add a product first!");
			    }
			else {
			secondFrameView s1=new secondFrameView(list);
			s1.setVisible(true);
			}
		}
	}
	
	class TreatButtonDeleteProduct implements ActionListener{
		public void actionPerformed(ActionEvent e) {
			if(list.isEmpty()==true) { 
				JOptionPane.showMessageDialog(null, 
						"There are no products! Please add a product first, then try to delete it!");
			    }
			else {
			secondFrameDelete s2=new secondFrameDelete(list);
			s2.setVisible(true);
			}
		}
	}
	
	class TreatButtonChangeQuantity implements ActionListener{
		public void actionPerformed(ActionEvent e) {
			if(list.isEmpty()==true) { 
				JOptionPane.showMessageDialog(null, 
				"There are no products! Please add a product first, then try to change the quantity!");
			    }
			else {
			secondFrameChangeQuantity s3=new secondFrameChangeQuantity(list);
			s3.setVisible(true);
			}
		}
	}
	
}

class Product{
	String name;
	String quantity;
	String price;
	public Product(String name, String quantity, String price) {
		this.name=name;
		this.quantity=quantity;
		this.price=price;
	}
	
	public Product(String name) {
		this.name=name;
	}
	
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getQuantity() {
		return quantity;
	}
	public void setQuantity(String quantity) {
		this.quantity = quantity;
	}
	public String getPrice() {
		return price;
	}
	public void setPrice(String price) {
		this.price = price;
	}
	@Override
	public String toString() {
		return "Product [name=" + name + ", quantity=" + quantity + ", price=" + price + "]";
	}
	
	
}

class secondFrameAdd extends JFrame{
	private JLabel name= new JLabel("Enter the name of the product:");
	private JLabel quantityAvailable= new JLabel("Enter the available quantity of the product:");
	private JLabel price= new JLabel("Enter the price of the product:");
    private JTextField productName= new JTextField(100);
    private JTextField productQuantity= new JTextField(100);
    private JTextField productPrice = new JTextField(100);
    private JButton addButton= new JButton("Add");
    
    public ArrayList<Product> list = new ArrayList<Product>();
	
	public secondFrameAdd(ArrayList<Product> list) {
		super("Add Product");
		this.list=list;
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setSize(450,400);
		initAdd();
		setVisible(true);
	}
	public void initAdd() {
        this.setLayout(null);
        name.setBounds(10, 10, 180, 100);
        add(name);
        productName.setBounds(200,50,200,20);
        add(productName);
        quantityAvailable.setBounds(10,80, 250, 100);
        add(quantityAvailable);
        productQuantity.setBounds(270,120,130, 20);
        add(productQuantity);
        price.setBounds(10, 200, 180, 20);
        add(price);
        productPrice.setBounds(200,200,200, 20);
        add(productPrice);      
        addButton.setBounds(160, 270, 100, 40);
        add(addButton);
        addButton.addActionListener(new TreatButtonAdd());
	}
	
	class TreatButtonAdd implements ActionListener{
		public void actionPerformed(ActionEvent e) {
			String sName,sQuantity,sPrice;
			int ok=1;
            sName = productName.getText();
			sQuantity=productQuantity.getText();
			sPrice=productPrice.getText();
	        boolean numericQuantity = true, numericPrice=true;
	        try {
	            Double num = Double.parseDouble(sQuantity);
	        } catch (NumberFormatException et) {
	            numericQuantity = false;
	        }
	        if(numericQuantity==false)
	        	JOptionPane.showMessageDialog(null, "Quantity must be a number! Only digits!");
	        try {
	            Double num1 = Double.parseDouble(sPrice);
	        } catch (NumberFormatException et) {
	            numericPrice = false;
	        }
	        if(numericPrice==false)
	        	JOptionPane.showMessageDialog(null, "Price must be a number! Only digits!");
	        
	        if(numericQuantity==true && numericPrice==true) {			
			for(int i=0; i<list.size();i++) {
				if(list.get(i).name.equals(sName)) {
					JOptionPane.showMessageDialog(null, "Product already exists!");
					ok=0;
				}
			}
			if(ok==1) {
			list.add(new Product(sName,sQuantity,sPrice));
			JOptionPane.showMessageDialog(null, "Product added!");
			}
	      }
		}
	}
}

class secondFrameView extends JFrame{
	private JTextArea textArea= new JTextArea();
	public ArrayList<Product> list=new ArrayList<Product>();
	
	public secondFrameView(ArrayList<Product> list) {
		super("View all products");
		this.list=list;
		setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		setSize(450,300);
		initView();
		setVisible(true);
	}
	public void initView() {
        this.setLayout(null);
        textArea.setBounds(10, 10, 410, 230);
        add(textArea);
        try {
      	  for(Product l:list) {
      	       secondFrameView.this.textArea.append(l.toString());
      	       secondFrameView.this.textArea.append("\n");
      	  }
        }catch(Exception ep){
  	          ep.printStackTrace();
      }
	}
}

class secondFrameDelete extends JFrame{
	public ArrayList<Product> list=new ArrayList<Product>();
	
	private JLabel name= new JLabel("Enter the name of the product:");
    private JTextField productName= new JTextField(100);
    private JButton deleteButton= new JButton("Delete");
    
	public secondFrameDelete(ArrayList<Product> list) {
		super("Delete product");
		this.list=list;		
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setSize(450,400);
			initDelete();
		    setVisible(true);		
	}
	
	public void initDelete() {
        this.setLayout(null);
        name.setBounds(10, 10, 180, 100);
        add(name);
        productName.setBounds(200,50,200,20);
        add(productName);
        deleteButton.setBounds(160, 270, 100, 40);
        add(deleteButton);
        deleteButton.addActionListener(new TreatButtonDelete());
	}
	
	class TreatButtonDelete implements ActionListener{
		public void actionPerformed(ActionEvent e) {
			String sName;
			int control=-1;
            sName = productName.getText();
			for(int i=0; i<list.size();i++) {
				if(list.get(i).name.equals(sName)){
					list.remove(i);
					control=1;
					JOptionPane.showMessageDialog(null, "Product deleted!");
				}
			}
			if(control==-1) JOptionPane.showMessageDialog(null, "Product was not found!");
		}
	}
}

class secondFrameChangeQuantity extends JFrame{
	public ArrayList<Product> list=new ArrayList<Product>();
	
	private JLabel name= new JLabel("Enter the name of the product:");
    private JTextField productName= new JTextField(100);
	private JLabel quantityAvailable= new JLabel("Enter the new quantity of the product:");
    private JTextField productQuantity= new JTextField(100);
    private JButton changeButton= new JButton("Change");
	
	public secondFrameChangeQuantity(ArrayList<Product> list) {
		super("Change quantity of a product");
		this.list=list;		
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			setSize(450,400);
			initChangeQuantity();
		    setVisible(true);		
	}
	
	public void initChangeQuantity() {
        this.setLayout(null);        
        name.setBounds(10, 10, 180, 100);
        add(name);
        productName.setBounds(200,50,200,20);
        add(productName);
        quantityAvailable.setBounds(10,80, 250, 100);
        add(quantityAvailable);
        productQuantity.setBounds(270,120,130, 20);
        add(productQuantity);
        changeButton.setBounds(160, 270, 100, 40);
        add(changeButton);
        changeButton.addActionListener(new TreatButtonChange());        
	}
	
	class TreatButtonChange implements ActionListener{
		public void actionPerformed(ActionEvent e) {
			String sName,sQuantity;
            sName = productName.getText();
            int control=-1;
            sQuantity=productQuantity.getText();
	        boolean numericQuantity = true;
	        try {
	            Double num = Double.parseDouble(sQuantity);
	        } catch (NumberFormatException et) {
	            numericQuantity = false;
	        }
	        if(numericQuantity==false)
	        	JOptionPane.showMessageDialog(null, "Quantity must be a number! Only digits!");
	        
	        if(numericQuantity==true) { 
			for(int i=0; i<list.size();i++) {
				if(list.get(i).name.equals(sName)) {
					list.get(i).quantity=sQuantity;
					control=1;
					JOptionPane.showMessageDialog(null, "Product quantity changed!");
				}
			}
			if(control==-1) JOptionPane.showMessageDialog(null, "Product was not found!");
	      }
		}
	}
}