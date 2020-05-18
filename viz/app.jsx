import React from "react";
import ReactDOM from "react-dom";
import Navbar from "react-bootstrap/Navbar";
import Nav from "react-bootstrap/Nav";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import { json } from "d3";
import PDP from "alignment.js/prevent_default_patch";

import Superreads from "./superreads.jsx";

PDP(document);

import "./style.scss";


function Header(props) {
  return (<div>
    <Navbar bg="light">
      <Navbar.Brand>SuperSEAL</Navbar.Brand>
    </Navbar>
    <Container>
      {props.children}
    </Container>
  </div>);
}

class App extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      superreads: null,
      viewing: "superreads"
    }
  }
  componentDidMount() {
    json("/superreads.json")
      .then(json => {
        json.number_of_sites = json.map(superread => superread.cv_end)
          .reduce((acc, curr) => Math.max(acc, curr), 0);
        this.setState({superreads: json});
      });
  }
  render() {
    switch (this.state.viewing) {
      case "superreads":
        return(<Header>
          <Superreads json={this.state.superreads} />
        </Header>);
    }
  }
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement("div"))
);
