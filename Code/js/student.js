const Student = {
  // please fill in your name and NetID
  // your NetID is the part of your email before @princeton.edu
  name: "Jyoti Sunkara",
  netID: "jyoti.sunkara@students.iiit.ac.in",
};

Student.updateHTML = function() {
  const studentInfo = this.name + " &lt;" + this.netID + "&gt;";
  document.getElementById("student").innerHTML = studentInfo;
};
