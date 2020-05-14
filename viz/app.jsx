import React from "react";
import ReactDOM from "react-dom";
import Navbar from "react-bootstrap/Navbar";
import Nav from "react-bootstrap/Nav";
import Container from "react-bootstrap/Container";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import { json } from "d3";
import PDP from "alignment.js/prevent_default_patch";

PDP(document);

import "./style.scss";


class App extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      superreads: null
    }
  }
  componentDidMount() {
    json("/superreads.json")
      .then(json => {
        console.log(json);
      });
  }
  render() {
    return (<div>
      <Navbar bg="light">
        <Navbar.Brand>SuperSEAL</Navbar.Brand>
      </Navbar>
      <Container>
      </Container>
    </div>);
  }
}

ReactDOM.render(
  <App />,
  document.body.appendChild(document.createElement("div"))
);
